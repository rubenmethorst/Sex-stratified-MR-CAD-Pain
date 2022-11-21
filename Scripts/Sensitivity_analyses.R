##             SCRIPT TO PERFORM MR STUDY: Sensitivity analyses            ##
##            Ruben Methorst, LUMC, 2021                      ##

VERSION="v1.0"
LASTEDITDATE="2022-11-21"
SCRIPTNAME="Sensitivity analyses for MR project, LUMC"
AUTHOR="Ruben Methorst | r.methorst@lumc.nl"
THISYEAR = format(as.Date(as.POSIXlt(Sys.time())), "%Y")

# Check if arguments are correct
args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("One argument must be supplied: directory were RDS resides", call.=FALSE)
} else {
  cat(paste0("\n\nOuput directory: ", args[1]))
}

# Set directory
path <- args[1]
path.to.data <- args[2]

# Data that should be in this directory
# - RDS Object from MR_analysis.R

# Data generated throughout this script will be saved in the same directory

# Opening statement
cat(paste0(
  "\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
",SCRIPTNAME,"
",VERSION," - ",LASTEDITDATE,"
",AUTHOR," | 1998-",THISYEAR,".
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))

# Lets time who is faster in analyzing the data
timeStart <- Sys.time()
date <- Sys.Date()

## Lets first install the required packages
# FUNCTION TO INSTALL PACKAGES: Courtsey of S.W. vanderlaan (swvanderlaan.github.io)
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented. 
    #update.install.packages.auto(ask = FALSE) 
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"https://cloud.r-project.org/\")", x)))
  }
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    BiocManager::install() # this would entail updating installed packages, which in turned may not be warrented
    
    eval(parse(text = sprintf("BiocManager::install(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

# FUNCTION TO DO THINGS QUIETLY
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Required packages
quiet(
  suppressPackageStartupMessages(c(
    install.packages.auto("tidyverse"),
    install.packages.auto("data.table"),
    install.packages.auto("cowplot"),
    install.packages.auto("crayon"),
    install.packages.auto("TwoSampleMR"),
    install.packages.auto("MRInstruments"),
    install.packages.auto("ieugwasr"),
    install.packages.auto("ggsignif"),
    install.packages.auto("readxl"),
    install.packages.auto("MRPRESSO")
  )))

cat(bold(green("Required packages are installed and loaded!\n\nLoading data...\n")))


## Load the data
dataset <- readRDS(file = paste0(path, "/MR_results_R_0.001_fstat_10_excl_LPA_outlier.rds"))

# min F-stat in all analyses
for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    cat("\n\nTesting for weak instrument bias\n")
    
    cat(paste0(names(dataset[i])), "-->",
        paste0(names(dataset[[i]])[j]),
        "f-stat:", min(dataset[[i]][[j]]$instruments_exp$f_stat),
        "\n")
    
    # F-stat < 20 is considered weak. We used a less stringet GWAS p-value, so more stringent on F-stat
    if (min(dataset[[i]][[j]]$instruments_exp$f_stat) < 10){
      cat("Weak instrument bias detected!!\n")
    }
  }
}


# significant MR Egger non-intercept
for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    try({
      if (dataset[[i]][[j]]$pleiotropy$pval < 0.05){
        cat(paste0(names(dataset[i])), "-->",
            paste0(names(dataset[[i]])[j]), 
            "significant pleiotropy!!\n")
      } else {
        cat(paste0(names(dataset[i])), "-->",
            paste0(names(dataset[[i]])[j]), 
            ": All Good!!\n")
        
      }
    })
  }
}


# Q-stat in all analyses
for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    try({
      # Perform het testing
      dataset[[i]][[j]]$heterogeneity <- mr_heterogeneity(dataset[[i]][[j]]$data)
      
      cat("\n\nTesting for weak instrument bias\n")
      
      if (dataset[[i]][[j]]$heterogeneity$Q_pval[2] < 0.05){
        cat("Heterogeneity detected!!\nFor:",
            paste0(names(dataset[i])), "-->",
            paste0(names(dataset[[i]])[j]))
      }
    })
  }
}

# Leave-one-out plots
# Make a single PDF file with all LOOs and titles
pdf(file = paste0(path, "LOO_plots.pdf"), onefile = T, width = 7.0, height = 7.5)

for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    try({
      # make a LOO
      res_loo <- mr_leaveoneout(dataset[[i]][[j]]$data)
      plot <-  mr_leaveoneout_plot(res_loo)
      
      print(plot)
    })
  }
}

dev.off()

# Directionality testing
for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    try({
      # Perform directionality test
      dataset[[i]][[j]]$directionality <- directionality_test(dataset[[i]][[j]]$data)
      
      if (dataset[[i]][[j]]$directionality$correct_causal_direction == "FALSE"){
        cat("Wrong direction!!\nFor:",
            paste0(names(dataset[i])), "-->",
            paste0(names(dataset[[i]])[j]))
        
      }
    })
  }
}


# Outlier analysis and re-running MR
# MR-PRESSO

for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    try({
      
      # Perform MR PRESSO
      dataset[[i]][[j]]$MR_PRESSO <- mr_presso(data = dataset[[i]][[j]]$data, 
                                               BetaOutcome = "beta.outcome", 
                                               BetaExposure = "beta.exposure", 
                                               SdOutcome = "se.outcome", 
                                               SdExposure = "se.exposure", 
                                               OUTLIERtest = TRUE, 
                                               DISTORTIONtest = TRUE, 
                                               NbDistribution = 3000, # STANDARD IS 1000!, but gave an accurarcy error
                                               SignifThreshold = 0.05)
      
      # Print the results
      cat(paste0(names(dataset[i])), "-->",
          paste0(names(dataset[[i]])[j]), "\n",
          dataset[[i]][[j]]$MR_PRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, "\t p-val hor. pleiotropy:",
          dataset[[i]][[j]]$MR_PRESSO$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
    })
  }
}



# Run analyses without MR-PRESSO identified outliers
for (i in 1:length(dataset)){
  for (j in 1:length(names(dataset[[i]]))){
    
    outliers <- c(dataset[[i]][[j]]$MR_PRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
    
    if (is.numeric(outliers) == TRUE){
      
      # Remove outliers and make new dataset
      exp_vs_out <- dataset[[i]][[j]]$data[-c(outliers),]
      dataset[[i]][[j]]$outlier_analysis$data <- exp_vs_out
      
      # Perform analyses without outliers
      dataset[[i]][[j]]$outlier_analysis$mr <- mr(exp_vs_out)
      dataset[[i]][[j]]$outlier_analysis$mr <- generate_odds_ratios(dataset[[i]][[j]]$outlier_analysis$mr)
      dataset[[i]][[j]]$outlier_analysis$single_snp <- mr_singlesnp(exp_vs_out)
      dataset[[i]][[j]]$outlier_analysis$pleiotropy <- mr_pleiotropy_test(exp_vs_out)
      
    }
  }
}


# PERFORM COMPARISON BETWEEN RESULTS (Interaction p-value)
# How to compare beta+SE: https://www.bmj.com/content/326/7382/219##
# Statistics to get p-value: https://www.bmj.com/content/343/bmj.d2304

for (i in 1:length(dataset[[1]])){
  
  diff <- dataset[[1]][[i]]$mr$b - dataset[[2]][[i]]$mr$b
  diff.se <- sqrt(dataset[[1]][[i]]$mr$se^2 + dataset[[2]][[i]]$mr$se^2)  
  z <- abs(diff / diff.se)
  p <- exp(-0.717*z - 0.417*z^2)
  
  cat(blue("Comparison male vs female: ",
           paste0(names(dataset[[2]])[i]), "\nP-value of interaction",
           p[3], "\t"))
  
  if (p[3] < 0.05){
    cat(bold(green("Significant")))
  }
  
  cat("\n\n")
  
  # Save the data
  dataset[[1]][[i]]$mr$interaction_sex_pval <- p
  dataset[[2]][[i]]$mr$interaction_sex_pval <- p
  
}

# Chest pain vs R07.4 comparison
# Male and Female CAD
for (i in c(1,2)){
  
  diff <- dataset[[i]][[1]]$mr$b - dataset[[i]][[5]]$mr$b
  diff.se <- sqrt(dataset[[i]][[1]]$mr$se^2 + dataset[[i]][[5]]$mr$se^2)  
  z <- abs(diff / diff.se)
  p <- exp(-0.717*z - 0.417*z^2)
  
  cat(blue("Comparison: ",
           names(dataset)[i], "Chest pain vs R07.4 \nP-value of interaction",
           p[3], "\t"))
  
  if (p[3] < 0.05){
    cat(bold(green("Significant")))
  }
  
  cat("\n\n")
  
  # Save the data
  dataset[[i]][[1]]$mr$interaction_chest_pval <- p
  
}



# Save dataset with validation
saveRDS(dataset, file = "MR_results_R_0.001_fstat_10_exlc_LPA_outlier_VALIDATED.rds")




# End of script
timeEnd <- Sys.time()

cat(paste0(bold(green("Finished in: ", difftime(timeEnd, timeStart, units='mins'), " Minutes.", 
                      "\nResults and Plots can be found here:", path, bold(silver("\n\n\t\tByebye\n\n\n\n\n"))))))

cat(paste0(yellow(italic(
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
  "The MIT License (MIT)\n",
  "Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n",
  "    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n",
  "Reference: http://opensource.org.\n",
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))))

## Ruben Methorst, 1998-2021  ##
