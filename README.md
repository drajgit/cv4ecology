Exaplanation of cv4ecology_code_example.R

Script Objective: To calculate and visualize univariate and multivariate outliers in 
morphometric measurements. Following traits had been measured: 
1.Body Length (BL) 
2.Right Forewing Length (RFW)
3.Left Forewing Length (LFW)
4.Right Hindwing Length (RHW)
5.Left Hindwing Length (LHW)

#### WORKFLOW ########
Function multivariate outlier executes the following
1. imports an ods file with image names and trait measurements. Each row stands for an individual
2. Clean data
3. Calculate average FW_Area and HW_Area
4. Calculate robust mean and standard deviation for BL, FW_Area and HW_Area. All the outlier 
values are flagged as 2 
5. Perform Regression between trait pairs. Calculate the robust mean and standard deviation of 
residuals
6. Flag the outliers as 3
7. Visualize the regression and plot data ellipse on it with different flag categories color coded
