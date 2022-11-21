#!/bin/bash
# Pipeline script to run the MR analyses in one go

DIR_PATH=/Users/rubenmethorst/Downloads/MR_test/
PATH_TO_DATA=/Users/rubenmethorst/Downloads/MR_test/Sex_stratified_data 

cd $DIR_PATH

echo "Making Manhattan plots..."
Rscript Scripts/Manhattan_plots.R $DIR_PATH $PATH_TO_DATA

echo "Performing MR analysis..."
Rscript Scripts/MR_analysis.R $DIR_PATH $PATH_TO_DATA

echo "Performing sensitivity analyses..."
Rscript Scripts/Sensitivity_analyses.R $DIR_PATH

echo "Visualizing the results..."
Rscript Scripts/Visualization_Figures.R $DIR_PATH

echo "DONE!"

# Script ended