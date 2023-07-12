##             SCRIPT TO PERFORM MR STUDY: MANHATTAN PLOTS            ##
##            Ruben Methorst, LUMC, 2021                      ##

VERSION="v1.0"
LASTEDITDATE="2022-11-21"
SCRIPTNAME="Drawing Manhatten plots from sex-stratified GWAS data, LUMC"
AUTHOR="Ruben Methorst | r.methorst@lumc.nl"
THISYEAR = format(as.Date(as.POSIXlt(Sys.time())), "%Y")

# Check if arguments are correct
args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("Two argument must be supplied: output directory + directory to GWAS data", call.=FALSE)
} else {
  cat(paste0("\n\nOuput directory: ", args[1], "\nData directory:", args[2]))
}

# Set directory
path <- args[1]
path.to.data <- args[2]

# For internal use only
path <- "~/Documents/Studie/PhD/MR project"
path.to.data <- "~/Documents/Studie/PhD/MR project/Sex-stratified data/Input/"


# Data that should be in this directory
# - BOLT_LLM GWAS results of performed GWASs

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
    install.packages.auto("GWASTools")
  )))

cat(bold(green("Required packages are installed and loaded!\n\nLoading data...\n")))

## Read UKBB data
# List files
files <- list.files(path = path.to.data, pattern = "*.txt", full.names = T)
names <- list.files(path = path.to.data, pattern = "*.txt", full.names = F)
names <- sub(path.to.data, "", names)

# Read the outcome data
datalist <- list()

for (i in 1:(length(files)-2)){
  temp <- quiet(read_delim(files[i], 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE, show_col_types = FALSE))
  
  cat("Lenght of:", names[i], "\t",length(temp$SNP), "\n")
  
  datalist[[i]] <- temp
  
  rm(temp)
  gc()
}

# Read the exposure data
for (i in (length(files)-1):length(files)){
  temp <- quiet(read_delim(files[i], 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE, show_col_types = FALSE))
  
  cat("Lenght of:", names[i], "\t",length(temp$SNP), "\n")
  
  datalist[[i]] <- temp
  
  rm(temp)
  gc()
}


names(datalist) <- c(names)

cat(bold(green("Data loaded, now making Manhattanplots (~15 min on M1 MacBook Pro (single core))...\n\n")))

# MyManhattan (https://github.com/alfonsosaera/myManhattan)

myManhattan <- function(df, graph.title = "", highlight = NULL, highlight.col = "green",
                        col = c("goldenrod1", "gold1"), even.facet = FALSE, chrom.lab = NULL,
                        suggestiveline = 1e-06, suggestivecolor = "grey",
                        genomewideline = 5e-08, genomewidecolor = "grey23",
                        font.size = 12, axis.size = 1, significance = NULL, report = FALSE,
                        inf.corr = 0.95, y.step = 5, point.size = 1.75){
  myMin <- min(df$P_BOLT_LMM[df$P_BOLT_LMM != 0]) * inf.corr
  df$P_BOLT_LMM[df$P_BOLT_LMM == 0] <- myMin
  require(ggplot2)
  require(stats)
  y.title <- expression(-log[10](italic(p)))
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      chrom.col <- c(chrom.col, col[1:dif])
    }
  }
  # y.max <- floor(max(-log10(df$P_BOLT_LMM))) + 1
  # if (y.max %% 2 != 0){
  #   y.max <- y.max + 1
  # }
  if (!is.null(chrom.lab)){
    if (length(unique(df$CHR)) != length(chrom.lab)){
      warning("Number of chrom.lab different of number of chromosomes in dataset, argument ignored.")
    } else {
      df$CHR <- factor(df$CHR, levels = unique(df$CHR), labels=chrom.lab)
    }
  }
  g <- ggplot(df) +
    geom_point(aes(BP, -log10(P_BOLT_LMM), colour = as.factor(CHR)), size = point.size)
  if (!is.null(significance)){
    if (is.numeric(significance)){
      genomewideline <- significance
      suggestiveline <- genomewideline / 0.005
    } else if (significance == "Bonferroni"){
      BFlevel <- 0.05 / length(df$SNP)
      cat("Bonferroni correction significance level:", BFlevel, "\n")
      genomewideline <- BFlevel
      suggestiveline <- BFlevel / 0.005
    } else if (significance == "FDR"){
      df$fdr <- p.adjust(df$P_BOLT_LMM, "fdr")
      genomewideline <- 0.05
      suggestiveline <- FALSE
      y.title <- expression(-log[10](italic(q)))
      g <- ggplot(df) +
        geom_point(aes(BP, -log10(fdr), colour = as.factor(CHR)), size = point.size)
      if (!is.null(highlight)) {
        if (is.numeric(highlight)){
          highlight <- as.character(df$SNP[df$P_BOLT_LMM < highlight])
        }
        if (any(!(highlight %in% df$SNP))){
          warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
        } else {
          g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                              aes(BP, -log10(fdr), group=SNP, colour=SNP),
                              color = highlight.col, size = point.size)
          highlight <- NULL
          # y.max <- floor(max(-log10(df$fdr))) + 1
          # if (y.max %% 2 != 0){
          #   y.max <- y.max + 1
          # }
        }
      }
    }
  }
  if (even.facet){
    g <- g + facet_grid(.~CHR, scale = "free_x", switch = "x")
  } else {
    g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  }
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, 30),
                       breaks = seq(from = 0, to = 30, by = y.step)) +
    scale_x_continuous() +
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          panel.background = element_blank(),
          axis.line.y = element_line(size = axis.size, color = "black"),
          axis.ticks.y = element_line(size = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "Chromosome", y = y.title)
  if (!is.null(highlight)) {
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P_BOLT_LMM < highlight])
    }
    if (any(!(highlight %in% df$SNP))){
      warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
    } else {
      g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                          aes(BP, -log10(P_BOLT_LMM), group=SNP, colour=SNP),
                          color = highlight.col, size = point.size)
    }
  }
  if (suggestiveline){
    g <- g + geom_hline(yintercept = -log10(suggestiveline), color = suggestivecolor)
  }
  if (genomewideline){
    g <- g + geom_hline(yintercept = -log10(genomewideline), color = genomewidecolor)
  }
  if (report){
    if (significance == "FDR"){
      rep <- df[df$fdr < 0.05, ]
    } else if (significance == "Bonferroni"){
      rep <- df[df$P_BOLT_LMM < BFlevel, ]
    } else if (is.numeric(significance)){
      rep <- df[df$P_BOLT_LMM < significance, ]
    } else {
      cat("using default significance level, 5e-8")
      rep <- df[df$P_BOLT_LMM < 5e-8, ]
    }
    print(rep)
  }
  return(g)
}

# ## SOMEHOW THIS DOESNT WORK IN A LOOP??? SO LETS DO IT STUPID FOR NOW
# ## MOST LIKELY A RAM PROBLEM (DATALIST == 35GB+)!!
# for (i in 1:17){
#   cat(bold(green("Manhattenplot:", names[i], "...\n\n")))
#   jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
#   myManhattan(datalist[[names[i]]], graph.title = names[i])
#   dev.off()
# }
# cat(bold(green("Done!\n\n")))

# 17 TIMES THE SAME
i <- 1
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 2
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 3
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 4
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 5
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 6
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 7
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 8
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 9
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 10
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

# CAD
i <- 11
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())

i <- 12
jpeg(file = paste0(path, "manhattan_plots/", names[i],"_Manhattan.jpeg"), width = 1520, height = 1080, units = "px")
myManhattan(datalist[[names[i]]], graph.title = names[i])
quiet(dev.off())


# Now make some QQ plots
cat(bold(green("Done!, now making some QQ plots...\n\n")))

i <- 1
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 2
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 3
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 4
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 5
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 6
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 7
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 8
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 9
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())


i <- 10
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 11
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())

i <- 12
jpeg(file = paste0(path, "/gwas_plots/", names[i],"_QQ_plot.jpeg"), width = 1080, height = 1080, units = "px")
qqPlot(datalist[[names[i]]][["P_BOLT_LMM"]])
quiet(dev.off())


# End of script
timeEnd <- Sys.time()

cat(paste0(bold(green("Finished in: ", difftime(timeEnd, timeStart, units='mins'), " Minutes.", 
                      "\nManhatten Plots can be found here:", path, bold(silver("\n\n\t\tByebye\n\n\n\n\n"))))))

cat(paste0(yellow(italic(
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
  "The MIT License (MIT)\n",
  "Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n",
  "    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n",
  "Reference: http://opensource.org.\n",
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))))

## Ruben Methorst, 1998-2021  ##
