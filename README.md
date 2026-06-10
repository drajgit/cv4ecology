## Overview

This script identifies and visualizes both **univariate** and **multivariate outliers** in morphometric measurements of moth specimens using robust statistical methods and regression-based analyses.

The following traits are analysed:

- Body Length (`BL`)
- Right Forewing Area (`RFW_Area`)
- Left Forewing Area (`LFW_Area`)
- Right Hindwing Area (`RHW_Area`)
- Left Hindwing Area (`LHW_Area`)

## WORKFLOW 
Function multivariate_outlier executes the following
1. Imports an ods file (structure like that of input_sheet_example) with image names and trait measurements. Each row stands for an individual
2. Clean data
3. Calculate average FW_Area and HW_Area
4. Calculate robust mean and standard deviation for BL, FW_Area and HW_Area. All the outlier 
values are flagged as 2 
5. Perform Regression between trait pairs. Calculate the robust mean and standard deviation of 
residuals
6. Flag the outliers as 3
7. Visualize the scatterplots and plot confidence ellipses (99%, 95% and 90%) on it with different flag categories color coded
