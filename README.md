## Exploring sex differences in pain manifestation of coronary artery disease: a Mendelian Randomization study 
R. Methorst1, M.R.M. Jongbloed1,2, R. Noordam3, M.C. DeRuiter1

1. Department of Anatomy and Embryology, Leiden University Medical Centre, Leiden, The Netherlands
2. Department of Cardiology, Leiden University Medical Centre, Leiden, The Netherlands
3. Department of Internal Medicine, Section of Gerontology and Geriatrics, Leiden University Medical Centre, Leiden, The Netherlands



# Abstract

**Importance:** Clinical manifestation of atherogenic cardiovascular disease is considered to differ between men and women. The anatomy and physiology of cardiac pain have not been properly characterized with respect to sex-related differences in humans. 

**Objective:** To investigate the relationship between coronary artery disease (CAD) risk and pain localizations in both sexes using Mendelian Randomization (MR). 

**Design:** We performed two-sample MR on summary-level data from sex-stratified genome-wide association studies of CAD as well as chest, neck and shoulder, back, and facial pain using data from UK Biobank. 

**Setting, Participants:** Sex-stratified Mendelian randomization analyses on UK Biobank cohort (N > 450,000).

**Exposure(s) and Main outcome:** UK Biobank participants having CAD and participants experiencing pain in different locations to assess the relationship between CAD and pain localization in a sex-stratified manner.

**Results:** We identified 32 and 19 instrumental variables associated with CAD for men and women, respectively. Genetically-influenced CAD was associated with increased risk of self-reported chest pain in men (OR: 1.27, CI: 1.21;1.33) and women (OR: 1.44, CI: 1.20;1.73) with similar results for clinical chest pain diagnoses (men OR: 1.22, CI: 1.17;1.26; women OR: 1.31, CI: 1.18;1.46). CAD was associated with back pain in women only (OR: 1.35, CI: 1.03;1.66) compared to men (OR: 0.99, CI: 0.93;1.07) (p-value for interaction: 0.030). Neck and shoulder pain showed a difference between men (OR: 0.97, CI: 0.90;1.03) and women (OR: 1.22, CI: 0.91;1.63) (p-value for interaction: 0.041). Sensitivity analysis did not indicate results were biased by directional pleiotropy. 

**Conclusions and Relevance:** We observed CAD drives different pain manifestations in men and women. Current efforts suggested cardiac pain is not majorly limited to the chest but is also referred to other regions in women at the population level. Future research is required to identify what causes these differences at multiple levels of research. 

# Usage

* `pipeline.sh`: Used to perform all analyses presented in the manuscript
* `Manhattan_plot.R`: Used to make manhattan plots of sex-stratified GWAS data
* `MR_analysis.R`: Used to perform MR analysis
* `Sensitivity_analyses.R`: Used to perform all presented sensitivity analyses in the manuscript
* `Visualization_Figures.R`: Used to make Figure 2 and Figure 3 of the manuscript

GWAS summary statistics are available upon request: r.methorst@lumc.nl
