## Using Mendelian randomization to explore sex differences in pain manifestation of coronary artery disease 
R. Methorst1, M.R.M. Jongbloed1,2, R. Noordam3, M.C. DeRuiter1

1. Department of Anatomy and Embryology, Leiden University Medical Centre, Leiden, The Netherlands
2. Department of Cardiology, Leiden University Medical Centre, Leiden, The Netherlands
3. Department of Geriatrics, Leiden University Medical Centre, Leiden, The Netherlands



# Abstract

Background and Aim: Clinical manifestation of atherogenic cardiovascular disease differs between men and women. Men often present with classical symptoms i.e., chest 
pain radiating to the left arm and jaw. Women show a wide range of symptoms considered atypical, including but not limited to back pain and nausea. Anatomical, 
functional, psychological, and moleculair pathways feeding cardiac pain have not been properly characterized with respect to sex-related differences in human at 
population level. 

Methods and Results: We performed MR on sex-stratified GWAS of CAD as well as chest, neck and shoulder, back, and facial pain using the UK Biobank cohort. We 
identified 32 and 19 strong instrumental variables associated with CAD for men and women, respectively. MR analysis revealed a causal link between CAD and 
self-reported chest pain for men (OR: 1.33, CI: 1.25 – 1.41) and women (OR: 1.44, CI: 1.20 – 1.73) as well for clinical chest pain (men OR: 1.22, CI: 1.17 — 1.26; 
women OR: 1.31, CI: 1.18 — 1.46). CAD displayed a causal association with back pain specific to women (OR: 1.35, CI: 1.03 — 1.66) compared to men (interaction 
p-value: 0.036). Manifestation of neck and shoulder pain revealed a trend towards a difference in causality between men and women (interaction p-value: 0.051).

Discussion: Our study showed strong causality of CAD with manifestation of chest pain in both men and women. However, no stronger effect in men was found in chest 
pain contrary to clinical observations suggesting pain perception of the chest at population level could be similar between the sexes. Back pain, as well as neck and 
shoulder pain, showed a more causal effect in women compared to men. Pain perception of the heart is not majorly limited to the chest but is also reflected to other 
regions in women at population level. Future research is required to identify what drives these differences at multiple levels. 

# Usage

* `pipeline.sh`: Used to perform all analyses presented in the manuscript
* `Manhattan_plot.R`: Used to make manhattan plots of sex-stratified GWAS data
* `MR_analysis.R`: Used to perform MR analysis
* `Sensitivity_analyses.R`: Used to perform all presented sensitivity analyses in the manuscript
* `Visualization_Figures.R`: Used to make Figure 2 and Figure 3 of the manuscript

GWAS summary statistics are available upon request: r.methorst@lumc.nl
