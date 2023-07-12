##             SCRIPT TO PERFORM MR STUDY: Visualization            ##
##            Ruben Methorst, LUMC, 2021                      ##

VERSION="v1.0"
LASTEDITDATE="2022-11-21"
SCRIPTNAME="Figure 2 and Figure 3 for MR project, LUMC"
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
    install.packages.auto("readxl")
  )))

cat(bold(green("Required packages are installed and loaded!\n\nLoading data...\n")))


## Load the data
dataset <- readRDS(file = paste0(path, "/MR_results_R_0.001_fstat_10_exlc_LPA_outlier_5e6_VALIDATED_2023-06-29.rds"))

# Data wrangling to get plot input
names <- c("male_chestpain", "female_chestpain", "male_R07.4", "female_R07.4")

exposure <- c(rep("CAD", 4))
outcome <- c(rep("Chest Pain", 2), rep("Clinical Chest Pain", 2))
Sex <- c("Male", "Female", "Male", "Female")

OR <- c(dataset[["result_male_CAD"]][["Chestpain_2335_UKB_MEN.txt"]][["mr"]]$or[3],
        dataset[["result_female_CAD"]][["Chestpain_2335_UKB_WOMEN.txt"]][["mr"]]$or[3],
        dataset[["result_male_CAD"]][["R074_UKB_MEN_European.txt"]][["mr"]]$or[3],
        dataset[["result_female_CAD"]][["R074_UKB_WOMEN_European.txt"]][["mr"]]$or[3])

LCI <- c(dataset[["result_male_CAD"]][["Chestpain_2335_UKB_MEN.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_female_CAD"]][["Chestpain_2335_UKB_WOMEN.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_male_CAD"]][["R074_UKB_MEN_European.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_female_CAD"]][["R074_UKB_WOMEN_European.txt"]][["mr"]]$or_lci95[3])

UCI <- c(dataset[["result_male_CAD"]][["Chestpain_2335_UKB_MEN.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_female_CAD"]][["Chestpain_2335_UKB_WOMEN.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_male_CAD"]][["R074_UKB_MEN_European.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_female_CAD"]][["R074_UKB_WOMEN_European.txt"]][["mr"]]$or_uci95[3])

Interaction <- c(dataset[["result_male_CAD"]][["Chestpain_2335_UKB_MEN.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_female_CAD"]][["Chestpain_2335_UKB_WOMEN.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_male_CAD"]][["R074_UKB_MEN_European.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_female_CAD"]][["R074_UKB_WOMEN_European.txt"]][["mr"]]$interaction_sex_pval[3])


# Extract chest pain data
chest_pain <- data.frame(names, exposure, outcome, Sex, OR, LCI, UCI, Interaction)

# Plot the Odds Ratio's
chest_pain_plot <- ggplot(chest_pain, aes(x = OR, 
                                          y = interaction(Sex, outcome))) +
  geom_point(aes(color=Sex), position = position_dodge(0.9), size = 4) + 
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = Sex, height = .25), size = 1.5) + 
  geom_vline(xintercept = 1, linetype="dotted") +
  geom_text(label="n.s.", x=1.05, y=1.5, color = "darkgray", size = 6) +
  geom_text(label="n.s.", x=1.05, y=3.5, color = "darkgray", size = 6) +
  geom_line(data = data.frame(x = c(1.1, 1.1), y = c(1,2)),
            aes(x = x, y = y), colour = "darkgray") +
  geom_line(data = data.frame(x = c(1.1, 1.1), y = c(3,4)),
            aes(x = x, y = y), colour = "darkgray") +
  xlab("Odds Ratio") + ggtitle("Sex-specific causal inference of chest pain manifestation in CAD") +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.spacing.x = unit(5,"mm"),
        legend.key.size = unit(2, 'cm'),
        legend.margin= margin(t = 5, unit = 'mm'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.margin = margin(t = 5, b = 5, r = 5, l = 5, unit='mm'),
        axis.title.x = element_text(size = 20, vjust = -1),
        axis.text.x = element_text(angle = 45, hjust = 1, size=14, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, size=16, colour = "black"),
        panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black"), aspect.ratio = 1) +
  scale_color_manual(values=c("salmon1", "khaki3", "salmon1", "khaki3"), labels = c("Women", "Men")) +
  scale_y_discrete(labels= c("Self-reported Chest Pain", "Self-reported Chest Pain", "Clinical Chest Pain", "Clinical Chest Pain")) +
  guides(colour = guide_legend(reverse=T))

chest_pain_plot

ggsave(filename = paste0(path, "/20230703_chest_pain_plot.pdf"),
       plot = chest_pain_plot + theme(plot.title = element_blank()),
       device = "pdf",
       width = 18*2, height = 15*2,
       units = "cm",
       dpi = 600)


# Chest Pain plot has been painted!
cat(bold(green("Chest Pain plot has been painted!...\n")))

# Data wrangling to get plot input
names <- c("male_back", "female_back", "male_neck", "female_neck", "male_facial", "female_facial")

exposure <- c(rep("CAD", 6))
outcome <- c(rep("Back Pain", 2), rep("Neck and Shoulder Pain", 2), rep("Facial Pain", 2))
Sex <- c("Male", "Female", "Male", "Female", "Male", "Female")

OR <- c(dataset[["result_male_CAD"]][["Pain_Back_UKB_MEN_European.txt"]][["mr"]]$or[3],
        dataset[["result_female_CAD"]][["Pain_Back_UKB_WOMEN_European.txt"]][["mr"]]$or[3],
        dataset[["result_male_CAD"]][["Pain_Neck_Schoulder_MEN_European.txt"]][["mr"]]$or[3],
        dataset[["result_female_CAD"]][["Pain_Neck_Schoulder_UKB_WOMEN_European.txt"]][["mr"]]$or[3],
        dataset[["result_male_CAD"]][["Pain_Facial_UKB_MEN_European.txt"]][["mr"]]$or[3],
        dataset[["result_female_CAD"]][["Pain_Facial_UKB_WOMEN_European.txt"]][["mr"]]$or[3])

LCI <- c(dataset[["result_male_CAD"]][["Pain_Back_UKB_MEN_European.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_female_CAD"]][["Pain_Back_UKB_WOMEN_European.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_male_CAD"]][["Pain_Neck_Schoulder_MEN_European.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_female_CAD"]][["Pain_Neck_Schoulder_UKB_WOMEN_European.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_male_CAD"]][["Pain_Facial_UKB_MEN_European.txt"]][["mr"]]$or_lci95[3],
         dataset[["result_female_CAD"]][["Pain_Facial_UKB_WOMEN_European.txt"]][["mr"]]$or_lci95[3])

UCI <- c(dataset[["result_male_CAD"]][["Pain_Back_UKB_MEN_European.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_female_CAD"]][["Pain_Back_UKB_WOMEN_European.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_male_CAD"]][["Pain_Neck_Schoulder_MEN_European.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_female_CAD"]][["Pain_Neck_Schoulder_UKB_WOMEN_European.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_male_CAD"]][["Pain_Facial_UKB_MEN_European.txt"]][["mr"]]$or_uci95[3],
         dataset[["result_female_CAD"]][["Pain_Facial_UKB_WOMEN_European.txt"]][["mr"]]$or_uci95[3])

Interaction <- c(dataset[["result_male_CAD"]][["Pain_Back_UKB_MEN_European.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_female_CAD"]][["Pain_Back_UKB_WOMEN_European.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_male_CAD"]][["Pain_Neck_Schoulder_MEN_European.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_female_CAD"]][["Pain_Neck_Schoulder_UKB_WOMEN_European.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_male_CAD"]][["Pain_Facial_UKB_MEN_European.txt"]][["mr"]]$interaction_sex_pval[3],
                 dataset[["result_female_CAD"]][["Pain_Facial_UKB_WOMEN_European.txt"]][["mr"]]$interaction_sex_pval[3])


# Extract chest pain data
other_pain <- data.frame(names, exposure, outcome, Sex, OR, LCI, UCI, Interaction)


# Plot the Odds Ratio's
other_pain_plot <- ggplot(other_pain, aes(x = OR, 
                                          y = interaction(Sex, outcome))) +
  geom_point(aes(color=Sex), position = position_dodge(0.9), size = 4) + 
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = Sex, height = .25), size = 1.5) + 
  geom_vline(xintercept = 1, linetype="dotted") +
  geom_text(label="*", x=0.67, y=1.5, color = "darkgray", size = 6) +
  geom_text(label="n.s.", x=0.67, y=3.5, color = "darkgray", size = 6) +
  geom_text(label="*", x=0.67, y=5.5, color = "darkgray", size = 6) +
  geom_line(data = data.frame(x = c(0.78, 0.78), y = c(1,2)),
            aes(x = x, y = y), colour = "darkgray") +
  geom_line(data = data.frame(x = c(0.78, 0.78), y = c(3,4)),
            aes(x = x, y = y), colour = "darkgray") +
  geom_line(data = data.frame(x = c(0.78, 0.78), y = c(5,6)),
            aes(x = x, y = y), colour = "darkgray") +
  xlab("Odds Ratio") + ggtitle("Sex-specific causal inference of non-typical pain manifestation in CAD") +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.spacing.x = unit(5,"mm"),
        legend.key.size = unit(2, 'cm'),
        legend.margin= margin(t = 5, unit = 'mm'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.margin = margin(t = 5, b = 5, r = 5, l = 5, unit='mm'),
        axis.title.x = element_text(size = 20, vjust = -1),
        axis.text.x = element_text(angle = 45, hjust = 1, size=14, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, size=16, colour = "black"),
        panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black"), aspect.ratio = 1) +
  scale_color_manual(values=c("salmon1", "khaki3", "salmon1", "khaki3","salmon1", "khaki3"), labels = c("Women", "Men")) +
  scale_y_discrete(labels= c(rep("Back Pain", 2), rep("Facial Pain", 2), rep("Neck and Shoulder Pain", 2))) +
  guides(colour = guide_legend(reverse=T)) + xlim(c(0.60, 2.0))

other_pain_plot

ggsave(filename = paste0(path, "/20230703_other_pain_plot.pdf"),
       plot = other_pain_plot + theme(plot.title = element_blank()),
       device = "pdf",
       width = 18*2, height = 17*2,
       units = "cm",
       dpi = 600)



# Other Pain plot has been painted!
cat(bold(green("Other Pain plot has been painted!...\n")))



# End of script
timeEnd <- Sys.time()

cat(paste0(bold(green("Finished in: ", difftime(timeEnd, timeStart, units='mins'), " Minutes.", 
                      "\nPlots can be found here:", path, bold(silver("\n\n\t\tByebye\n\n\n\n\n"))))))

cat(paste0(yellow(italic(
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
  "The MIT License (MIT)\n",
  "Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n",
  "    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n",
  "Reference: http://opensource.org.\n",
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))))

## Ruben Methorst, 1998-2021  ##
