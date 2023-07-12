##             SCRIPT TO PERFORM MR STUDY: MR Analysis            ##
##            Ruben Methorst, LUMC, 2021                      ##

VERSION="v1.0"
LASTEDITDATE="2022-11-21"
SCRIPTNAME="Perform sex-stratified MR CAD -> Pain localisation, LUMC"
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
# FUNCTION TO INSTALL PACKAGES: Courtsey of S.W. van der Laan (swvanderlaan.github.io)
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
    install.packages.auto("readxl")
  )))

cat(bold(green("Required packages are installed and loaded!\n\nLoading data...\n")))

## SET SOME IMPORTANT VARIABLES
# fstat:  cut-off f-statistic
# outlier_m/f: the SNP in the LPA locus introducing bias in MR
fstat <- 10
outlier_m <- "rs55730499"
outlier_f <- "rs118039278"


## Read UKBB data
# List files
files <- list.files(path = path.to.data, pattern = "*.txt", full.names = T)
names <- list.files(path = path.to.data, pattern = "*.txt", full.names = F)
names <- sub(path.to.data, "", names)

# Sample sizes
Sample_size <- read_excel(paste0(path, "/GitHub/Sex-stratified-MR-CAD-Pain/Sample_size.xlsx"))

# Read the outcome data
datalist <- list()

for (i in 1:(length(files)-2)){
  temp <- quiet(read_delim(files[i], 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE, show_col_types = FALSE))
  
  # Add sample sizes
  temp$case <- Sample_size$sample_size[[i]]
  temp$control <- Sample_size$controls[[i]]
  
  cat("Lenght of:", names[i], "\t",length(temp$SNP), "\n")
  
  datalist[[i]] <- temp
  
  rm(temp)
  quiet(gc())
}


## DO CAD seperately, female and male have a different GWAS cut-off (see methods paper)

# Read male CAD data
CAD_m <- quiet(read_delim(files[11], 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE, show_col_types = FALSE))

# Filter data
CAD_m <- subset(CAD_m, P_BOLT_LMM < 5e-6) 

# Add sample sizes
CAD_m$case <- Sample_size$sample_size[[11]]
CAD_m$control <- Sample_size$controls[[11]]

cat("Lenght of:", names[11], "\t",length(CAD_m$SNP), "\n")

datalist[[11]] <- CAD_m

rm(CAD_m)
quiet(gc())

# Read female CAD data
CAD_f <- quiet(read_delim(files[12], 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE, show_col_types = FALSE))

# Filter data
CAD_f <- subset(CAD_f, P_BOLT_LMM < 5e-6)

# Add sample sizes
CAD_f$case <- Sample_size$sample_size[[12]]
CAD_f$control <- Sample_size$controls[[12]]

cat("Lenght of:", names[12], "\t",length(CAD_f$SNP), "\n")

datalist[[12]] <- CAD_f

rm(CAD_f)
quiet(gc())


# Name de list of dataframes
names(datalist) <- c(names)

# Remove empty dataframes for less problems down the line
for (i in 1:length(datalist)){
  if (length(datalist[[i]]$SNP) == 0){
    datalist[[i]] <- NULL
  }
}

# Segregate in male and female datasets
males <- c("_M_", "_MEN")
females <- c("_W_", "_WOMEN")

male_data <- datalist[grepl(paste(males, collapse = "|"), names(datalist))]
female_data <- datalist[grepl(paste(females, collapse = "|"), names(datalist))]

rm(datalist)
quiet(gc())

# NOTE: Last two elements are CAD and MI

cat(bold(green("Data loaded, now onto MR...\n\n")))

# Perform MR
# 4 exposures and 2 outcomes

### MALES -------
## CAD as outcome
#read exposure data as MR input
write.csv(as.data.frame(male_data[[6]]), file = paste0(path, "/exposure_input.csv"))

exposure_data <- read_exposure_data(
  filename = paste0(path, "/exposure_input.csv"),
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM", 
  ncase_col = "case",
  ncontrol_col = "control"
)

exposure_data$exposure <- names(male_data)[6]

# Perform clumping
exposure_data <- clump_data(exposure_data, clump_r2 = 0.001)

# # Perform local clumping (NOT WORKING FOR NOW)
# exposure_data <- ld_clump(dplyr::tibble(rsid=exposure_data$SNP, pval=exposure_data$pval.exposure, id=exposure_data$exposure),
#                           clump_r2 = 0.001,
#                           plink_bin = genetics.binaRies::get_plink_binary(), 
#                           bfile = "~/Downloads/1kg.v3/EUR")

# Filter for F-stat
# # Standard F-statsicis
# F   = B^2/seBeta^2

# # F-stat according to Pierce et al. IJE 2011
# # R^2: Proportion of variability
# # n:   Sample size (22323 cases)
# # k:   Number of IVs (33 in the end)
B <- abs(exposure_data$beta.exposure)
seBeta <- exposure_data$se.exposure
k <- length(exposure_data$SNP)
n <- exposure_data$samplesize.exposure
rsqrd <- get_r_from_bsen(B, seBeta, n)

F = (rsqrd*(n - 1 - k)) / ((1 - rsqrd)*k)

# Add to data and filter for F-stat
exposure_data$f_stat <- F
exposure_data <- subset(exposure_data, exposure_data$f_stat > fstat)

# Remove LPA outlier
exposure_data <- subset(exposure_data, exposure_data$SNP != outlier_m)


# Loop for MR analyses
result_male_CAD <- list()

for (i in 1:5){
  try({
    # Outcome data in format for MR
    write_csv(as.data.frame(male_data[[i]]), file = paste0(path, "/outcome_input.csv"))
    
    outcome_dat <- read_outcome_data(
      snps = exposure_data$SNP,
      filename = paste0(path, "/outcome_input.csv"),
      sep = ",",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "P_BOLT_LMM", 
      ncase_col = "case",
      ncontrol_col = "control"
    )
    
    outcome_dat$outcome <- names(male_data)[i]
    
    #harmonize data
    
    exp_vs_out <- harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = outcome_dat, 
      action = 2 
    )
    
    # Steiger filter test
    steiger <- steiger_filtering(exp_vs_out)
    exp_vs_out <- subset(steiger, steiger_pval < 0.05)
    
    # MR
    res_out <- mr(exp_vs_out)
    res_single_out <- mr_singlesnp(exp_vs_out)
    res_pleiotropy <- mr_pleiotropy_test(exp_vs_out)
    
    print(res_single_out)
    
    #save to a list
    result_male_CAD[[names(male_data)[i]]][["mr"]] <- generate_odds_ratios(res_out)
    result_male_CAD[[names(male_data)[i]]][["singlesnp"]] <- res_single_out
    result_male_CAD[[names(male_data)[i]]][["pleiotropy"]] <- res_pleiotropy
    result_male_CAD[[names(male_data)[i]]][["data"]] <- exp_vs_out
    result_male_CAD[[names(male_data)[i]]][["instruments_exp"]] <- exposure_data
    result_male_CAD[[names(male_data)[i]]][["instruments_out"]] <- outcome_dat
    
  })
}

### FEMALES
## CAD as outcome
#read exposure data as MR input
write.csv(as.data.frame(female_data[[6]]), file = paste0(path, "/exposure_input.csv"))

exposure_data <- read_exposure_data(
  filename = paste0(path, "/exposure_input.csv"),
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM", 
  ncase_col = "case",
  ncontrol_col = "control"
)

exposure_data$exposure <- names(female_data)[6]

# Perform clumping
exposure_data <- clump_data(exposure_data, clump_r2 = 0.001)

# # Perform local clumping (NOT WORKING FOR NOW)
# exposure_data <- ld_clump(dplyr::tibble(rsid=exposure_data$SNP, pval=exposure_data$pval.exposure, id=exposure_data$exposure),
#                           clump_r2 = 0.001,
#                           plink_bin = genetics.binaRies::get_plink_binary(), 
#                           bfile = "~/Downloads/1kg.v3/EUR")

# Filter for F-stat
# # Standard F-statsicis
# F   = B^2/seBeta^2

# # F-stat according to Pierce et al. IJE 2011
# # R^2: Proportion of variability
# # n:   Sample size (22323 cases)
# # k:   Number of IVs (33 in the end)
B <- abs(exposure_data$beta.exposure)
seBeta <- exposure_data$se.exposure
k <- length(exposure_data$SNP)
n <- exposure_data$samplesize.exposure
rsqrd <- get_r_from_bsen(B, seBeta, n)

F = (rsqrd*(n - 1 - k)) / ((1 - rsqrd)*k)

# Add to data and filter for F-stat
exposure_data$f_stat <- F
exposure_data <- subset(exposure_data, exposure_data$f_stat > fstat)

# Remove LPA outlier
exposure_data <- subset(exposure_data, exposure_data$SNP != outlier_f)

# Loop for MR analyses
result_female_CAD <- list()

for (i in 1:5){
  try({
    # Outcome data in format for MR
    write_csv(as.data.frame(female_data[[i]]), file = paste0(path, "/outcome_input.csv"))
    
    outcome_dat <- read_outcome_data(
      snps = exposure_data$SNP,
      filename = paste0(path, "/outcome_input.csv"),
      sep = ",",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "P_BOLT_LMM", 
      ncase_col = "case",
      ncontrol_col = "control"
    )
    
    outcome_dat$outcome <- names(female_data)[i]
    
    #harmonize data
    
    exp_vs_out <- harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = outcome_dat
    )
    
    # Steiger filter test
    steiger <- steiger_filtering(exp_vs_out)
    exp_vs_out <- subset(steiger, steiger_pval < 0.05)
    
    # MR
    res_out <- mr(exp_vs_out)
    res_single_out <- mr_singlesnp(exp_vs_out)
    res_pleiotropy <- mr_pleiotropy_test(exp_vs_out)
    
    print(res_single_out)
    
    #save to a list
    result_female_CAD[[names(female_data)[i]]][["mr"]] <- generate_odds_ratios(res_out)
    result_female_CAD[[names(female_data)[i]]][["singlesnp"]] <- res_single_out
    result_female_CAD[[names(female_data)[i]]][["pleiotropy"]] <- res_pleiotropy
    result_female_CAD[[names(female_data)[i]]][["data"]] <- exp_vs_out
    result_female_CAD[[names(female_data)[i]]][["instruments_exp"]] <- exposure_data
    result_female_CAD[[names(female_data)[i]]][["instruments_out"]] <- outcome_dat
    
  })
}


# Make a list and save
dataset <- list(result_male_CAD, result_female_CAD)
names(dataset) <- c("result_male_CAD", "result_female_CAD")

saveRDS(dataset, file = paste0(path, "/MR_results_R_0.001_fstat_10_excl_LPA_outlier_", date, ".rds"))



# End of script
timeEnd <- Sys.time()

cat(paste0(bold(green("Finished in: ", difftime(timeEnd, timeStart, units='mins'), " Minutes.", 
                      "\nResults can be found here:", path, bold(silver("\n\n\t\tByebye\n\n\n\n\n"))))))

cat(paste0(yellow(italic(
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
  "The MIT License (MIT)\n",
  "Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n",
  "    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n",
  "Reference: http://opensource.org.\n",
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))))

## Ruben Methorst, 1998-2021  ##