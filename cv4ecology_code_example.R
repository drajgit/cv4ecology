# Script Objective: To calculate and visualize univariate and multivariate outliers in 
# morphometric measurements. Following traits had been measured: 
# 1.Body Length (BL) 
# 2.Right Forewing Length (RFW)
# 3.Left Forewing Length (LFW)
# 4.Right Hindwing Length (RHW)
# 5.Left Hindwing Length (LHW)

#### WORKFLOW ########
# Function multivariate outlier executes the following
# 1. imports an ods file with image names and trait measurements. Each row stands for an individual
# 2. Clean data
# 3. Calculate average FW_Area and HW_Area
# 4. Calculate robust mean and standard deviation for BL, FW_Area and HW_Area. All the outlier 
# values are flagged as 2 
# 5. Perform Regression between trait pairs. Calculate the robust mean and standard deviation of 
# residuals
# 6. Flag the outliers as 3
# 7. Visualize the regression and plot data ellipse on it with different flag categories color coded

library(readODS)
library(tidyr)
library(ggplot2)
library(dplyr)
library(astro)
library(RColorBrewer)
library(RanglaPunjab)
library(scales)
library(patchwork)

multivariate_outlier = function(file){
  fl = read_ods(file)
  # convert all character NAs to NA
  fl[fl == "NA"] = NA
  # remove non_moth entries
  fl = fl %>%
    filter(!grepl("non_moth", Family))
  # spreadsheet has column Image_Pair. Upperside and Underside image file names of an individual are 
  # separated by ";"
  fl = fl %>%
    separate(Image_Pair, into = c("Upperside", "Underside"), sep = ";", remove = FALSE)
  
  # convert trait values to numeric
  fl$BL = as.numeric(fl$BL)
  fl$RFW_Area = as.numeric(fl$RFW_Area)
  fl$LFW_Area = as.numeric(fl$LFW_Area)
  fl$RHW_Area = as.numeric(fl$RHW_Area)
  fl$LHW_Area = as.numeric(fl$LHW_Area)
  
  # input spreadsheet had a column with flag for each trait. Value 0 stands for measured and 1 for
  # not measured
  fl$BL_Flag = as.numeric(fl$BL_Flag)
  fl$RFW_Flag = as.numeric(fl$RFW_Flag)
  fl$LFW_Flag = as.numeric(fl$LFW_Flag)
  fl$RHW_Flag = as.numeric(fl$RHW_Flag)
  fl$LHW_Flag = as.numeric(fl$LHW_Flag)
  
  # create two columns: FW_Area and HW_Area
  fl$FW_Area = NA
  fl$HW_Area = NA
  
  # calculate total FW Area. Forewing was digitized in upperside/dorsal images. 
  for(i in 1:nrow(fl)){
    if(!is.na(fl$Upperside[i])){
      # if both forewings measured
      if(!is.na(fl$RFW_Flag[i]) && !is.na(fl$LFW_Flag[i]) &&
         fl$RFW_Flag[i] == 0 && fl$LFW_Flag[i] == 0){
        fl$FW_Area[i] = (fl$RFW_Area[i] + fl$LFW_Area[i]) / 2
      }
      # only right forewing measured
      if(!is.na(fl$RFW_Flag[i]) && fl$RFW_Flag[i] == 0 &&
         (is.na(fl$LFW_Flag[i]) || fl$LFW_Flag[i] == 1)){
        fl$FW_Area[i] = fl$RFW_Area[i]
      }
      # only left forewing measured
      if(!is.na(fl$LFW_Flag[i]) && fl$LFW_Flag[i] == 0 &&
         (is.na(fl$RFW_Flag[i]) || fl$RFW_Flag[i] == 1)){
        fl$FW_Area[i] = fl$LFW_Area[i]
      }
      
    }
  }
  
  # calculate total HW Area. Hindwings were digitized in underside/ventral images
  for(j in 1:nrow(fl)){
    if(!is.na(fl$Underside[j])){
      # both hindwings measured
      if(!is.na(fl$RHW_Flag[j]) && !is.na(fl$LHW_Flag[j]) &&
         fl$RHW_Flag[j] == 0 && fl$LHW_Flag[j] == 0){
        fl$HW_Area[j] = (fl$RHW_Area[j] + fl$LHW_Area[j]) / 2
      }
      # only right hindwing measured
      if(!is.na(fl$RHW_Flag[j]) && fl$RHW_Flag[j] == 0 &&
         (is.na(fl$LHW_Flag[j]) || fl$LHW_Flag[j] == 1)){
        fl$HW_Area[j] = fl$RHW_Area[j]
      }
      # only left hindwing measured
      if(!is.na(fl$LHW_Flag[j]) && fl$LHW_Flag[j] == 0 &&
         (is.na(fl$RHW_Flag[j]) || fl$RHW_Flag[j] == 1)){
        fl$HW_Area[j] = fl$LHW_Area[j]
      }
      
    }
  }
  # create a copy of fl to be used for further steps
  df = fl
  # new column that mentions row number
  df$row_id = 1:nrow(df)
  # convert the traits into logarithms. This removes the skewness and makes the trait distributions
  # more symmetric
  
  df$BL = log10(df$BL)
  
  df$FW_Area = log10(df$FW_Area)
  df$HW_Area = log10(df$HW_Area)
  
  # function to calculate robust mean and sd
  robust_stat = function(x, n){
    x = x[!is.na(x)]
    # remove non-finite values
    x = x[is.finite(x)]
    x_stat = scmean(x, mult = n, loop = 10)
    mean_x = x_stat$m
    sd_x = x_stat$s
    return(list(mean_x = mean_x,
                sd_x = sd_x))
    
  }
  
  stat_bl = robust_stat(df$BL, 3.5)
  stat_fw = robust_stat(df$FW_Area, 3.5)
  stat_hw = robust_stat(df$HW_Area, 3.5)
  
  # Flag everything outside 3.5 sigma as 2 in BL_Flag
  for(ii in 1:nrow(df)){
    
    if(!is.na(df$BL_Flag[ii]) && df$BL_Flag[ii] == 0 && !is.na(df$BL[ii])){
      
      if((df$BL[ii] <= stat_bl$mean_x - 3.5 * stat_bl$sd_x) | 
         (df$BL[ii] >= stat_bl$mean_x + 3.5 * stat_bl$sd_x)){
        
        df$BL_Flag[ii] = 2
        
      }
      
    }
    
  }
  
  
  # Create a new flag with FW outliers
  df$FW_Flag = NA
  for(jj in 1:nrow(df)){
    if(!is.na(df$FW_Area[jj])){
      df$FW_Flag[jj] = 0
    }
  }
  
  # Flag outliers for FW
  for(f in 1:nrow(df)){
    if(!is.na(df$FW_Flag[f]) && df$FW_Flag[f] == 0 && !is.na(df$FW_Area[f])){
      
      if((df$FW_Area[f] <= stat_fw$mean_x - 3.5 * stat_fw$sd_x) | 
         (df$FW_Area[f] >= stat_fw$mean_x + 3.5 * stat_fw$sd_x)){
        
        df$FW_Flag[f] = 2
        
      }
      
    }
  }
  
  # Create a new flag with HW outliers
  df$HW_Flag = NA
  for(h in 1:nrow(df)){
    if(!is.na(df$HW_Area[h])){
      df$HW_Flag[h] = 0
    }
  }
  
  # create outlier flags for HW
  for(hh in 1:nrow(df)){
    if(!is.na(df$HW_Flag[hh]) && df$HW_Flag[hh] == 0 && !is.na(df$HW_Area[hh])){
      
      if((df$HW_Area[hh] <= stat_hw$mean_x - 3.5 * stat_hw$sd_x) | 
         (df$HW_Area[hh] >= stat_hw$mean_x + 3.5 * stat_hw$sd_x)){
        
        df$HW_Flag[hh] = 2
        
      }
      
    }
    
  }
  
  # REGRESSION BETWEEN BL and FW_AREA
  
  df01 <- df %>%
    filter(!is.na(BL_Flag) & BL_Flag %in% c(0,2))
  
  df02 = df01 %>%
    filter(!is.na(FW_Flag) & FW_Flag %in% c(0,2))
  
  model1 = lm(FW_Area~BL, data = df02)
  residuals1 = resid(model1)
  
  # calculate robust mean and sd of residuals
  res_stat1 = scmean(residuals1)
  res_mean1 = res_stat1$m
  res_sigma1 = res_stat1$s
  
  lower_cut1 = res_mean1 - 3.5 * res_sigma1
  upper_cut1 = res_mean1 + 3.5 * res_sigma1
  
  # Mark rows as  where the residual is too small OR too large
  outlier_rows = residuals1 < lower_cut1 | residuals1 > upper_cut1
  
  # rows where outlier_rows is true are selected and flag value for the two traits are updated to 3
  df02$FW_Flag[outlier_rows] = 3
  df02$BL_Flag[outlier_rows] = 3
  
  # Update the flag value for BL and FW in original df
  df$BL_Flag[df$row_id %in% df02$row_id] =
    pmax(df$BL_Flag[df$row_id %in% df02$row_id], df02$BL_Flag)
  df$FW_Flag[df$row_id %in% df02$row_id] =
    pmax(df$FW_Flag[df$row_id %in% df02$row_id], df02$FW_Flag)
  
  # Create a new column in df02 to assign legend outlier kind values
  df02$Outlier = "Good"
  df02$Outlier[df02$FW_Flag == 2] = "FW_Area_Outlier"
  df02$Outlier[df02$BL_Flag == 2] = "BL_Outlier"
  df02$Outlier[df02$FW_Flag == 3] = "Regression_Outlier"
  df02$Outlier[df02$BL_Flag == 3] = "Regression_Outlier"
  df02$Outlier[df02$FW_Flag == 2 & df02$BL_Flag == 2] = "BL_FW_Area_Outlier" 
  
  # create a vector of colors for visualization
  rangla_panjab <- c(
    "#0072B2",  # deep blue
    "#56B4E9",  # sky blue
    "#009E73",  # green
    "#A93226",  # yellow
    "#E69F00",  # mustard
    "#D55E00",  # vermillion
    "#CC79A7",  # magenta
    "#7F7F7F",  # grey
    "#332288",  # indigo
    "#88CCEE",  # light blue
    "#B10DC9"   # bright pink
  )
  rangla_panjab2 = RanglaPunjab("Phulkari")
  outlier_colors <- c(
    "Good" = "#0072B2",
    "BL_Outlier" = "red4",
    "FW_Area_Outlier" = rangla_panjab[3],
    "HW_Area_Outlier" = rangla_panjab[3],
    "BL_FW_Area_Outlier" = "black",
    "BL_HW_Area_Outlier" = rangla_panjab2[4],
    "FW_HW_Area_Outlier" = rangla_panjab2[4],
    "Regression_Outlier" = rangla_panjab[11]
  )
  p1 = ggplot(df02, aes(x = BL , y = FW_Area)) +
    geom_point(aes(color = Outlier),
               size = 0.6, alpha = 0.6) +
    
    stat_ellipse(aes(group = 1, linetype = "99%"),
                 type = "norm",
                 level = 0.99,
                 color = rangla_panjab[9],
                 linewidth = 0.8) +
    stat_ellipse(aes(group = 1, linetype = "95%"),
                 type = "norm",
                 level = 0.95,
                 color = "deeppink4",
                 linewidth = 0.8) +
    stat_ellipse(aes(group = 1, linetype = "90%"),
                 type = "norm",
                 level = 0.90,
                 color = rangla_panjab[5],
                 linewidth = 0.8) +
    scale_color_manual(
      name = "Outlier Category",
      values = outlier_colors,
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = "Confidence ellipse",
      values = c(
        "90%" = "solid",
        "95%" = "dashed",
        "99%" = "dashed"
      )
    )+
    guides(
      color = guide_legend(
        title = "Outlier Category",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 2)
      ),
      linetype = guide_legend(
        title = "Confidence ellipse",
        override.aes = list(linewidth = 1.2)
      )
    )+
    xlab(expression(log[10]("Body length (mm)"))) +
    ylab(expression(log[10](Forewing~area~(mm^2)))) +
    coord_fixed(xlim = c(0,3), ylim = c(0,3)) +
    
    theme_bw()+
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.key.size = unit(0.35, "cm")
    )
  
  #################################################################################################
  
  # REGRESIION BETWEEN BL AND HW_AREA
  df03 <- df %>%
    filter(!is.na(BL_Flag) & BL_Flag %in% c(0,2))
  
  df04 = df03 %>%
    filter(!is.na(HW_Flag) & HW_Flag %in% c(0,2))
  
  model2 = lm(HW_Area~BL, data = df04)
  residuals2 = resid(model2)
  
  # calculate robust mean and sd of residuals
  res_stat2 = scmean(residuals2)
  res_mean2 = res_stat2$m
  res_sigma2 = res_stat2$s
  
  lower_cut2 = res_mean2 - 3.5 * res_sigma2
  upper_cut2 = res_mean2 + 3.5 * res_sigma2
  
  # Mark rows as  where the residual is too small OR too large
  outlier_rows = residuals2 < lower_cut2 | residuals2 > upper_cut2
  
  # rows where outlier_rows is true are selected and flag value for the two traits are updated to 3
  df04$HW_Flag[outlier_rows] = 3
  df04$BL_Flag[outlier_rows] = 3
  
  # Update the flag value for BL and FW in original df
  df$BL_Flag[df$row_id %in% df04$row_id] =
    pmax(df$BL_Flag[df$row_id %in% df04$row_id], df04$BL_Flag)
  df$HW_Flag[df$row_id %in% df04$row_id] =
    pmax(df$HW_Flag[df$row_id %in% df04$row_id], df04$HW_Flag)
  
  # Cretae a new column in df02 to assign legend outlier kind values
  df04$Outlier = "Good"
  df04$Outlier[df04$HW_Flag == 2] = "HW_Area_Outlier"
  df04$Outlier[df04$BL_Flag == 2] = "BL_Outlier"
  df04$Outlier[df04$HW_Flag == 3] = "Regression_Outlier"
  df04$Outlier[df04$BL_Flag == 3] = "Regression_Outlier"
  df04$Outlier[df04$HW_Flag == 2 & df04$BL_Flag == 2] = "BL_HW_Area_Outlier" 
  
  p2 = ggplot(df04, aes(x = BL , y = HW_Area)) +
    geom_point(aes(color = Outlier),
               size = 0.6, alpha = 0.6) +
    stat_ellipse(aes(group = 1, linetype = "99%"),
                 type = "norm",
                 level = 0.99,
                 color = rangla_panjab[9],
                 linewidth = 0.8) +
    stat_ellipse(aes(group = 1, linetype = "95%"),
                 type = "norm",
                 level = 0.95,
                 color = "deeppink4",
                 linewidth = 0.8) +
    stat_ellipse(aes(group = 1, linetype = "90%"),
                 type = "norm",
                 level = 0.90,
                 color = rangla_panjab[5],
                 linewidth = 0.8) +
    scale_color_manual(
      name = "Outlier Category",
      values = outlier_colors,
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = "Confidence ellipse",
      values = c(
        "90%" = "solid",
        "95%" = "dashed",
        "99%" = "dashed"
      )
    )+
    guides(
      color = guide_legend(
        title = "Outlier Category",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 2)
      ),
      linetype = guide_legend(
        title = "Confidence ellipse",
        override.aes = list(linewidth = 1.2)
      )
    )+
    xlab(expression(log[10]("Body length (mm)"))) +
    ylab(expression(log[10](Hindwing~area~(mm^2)))) +
    
    coord_fixed(xlim = c(0,3), ylim = c(0,3))+
    theme_bw()+
    theme(
      legend.position = "none",   # ← ADD THIS
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.key.size = unit(0.35, "cm")
    )
  
  
  
  ##################################################################################################
  # REGRESSION BETWEEN FW_AREA AND HW_AREA
  df05 <- df %>%
    filter(!is.na(FW_Flag) & FW_Flag %in% c(0,2))
  
  df06 = df05 %>%
    filter(!is.na(HW_Flag) & HW_Flag %in% c(0,2))
  
  model3 = lm(HW_Area~FW_Area, data = df06)
  residuals3 = resid(model3)
  
  res_stat3 = scmean(residuals3)
  res_mean3 = res_stat3$m
  res_sigma3 = res_stat3$s
  
  lower_cut3 = res_mean3 - 3.5 * res_sigma3
  upper_cut3 = res_mean3 + 3.5 * res_sigma3
  
  # Mark rows as  where the residual is too small OR too large
  outlier_rows = residuals3 < lower_cut3 | residuals3 > upper_cut3
  
  df06$HW_Flag[outlier_rows] = 3
  df06$FW_Flag[outlier_rows] = 3
  
  df$HW_Flag[df$row_id %in% df06$row_id] =
    pmax(df$HW_Flag[df$row_id %in% df06$row_id], df06$HW_Flag)
  df$FW_Flag[df$row_id %in% df06$row_id] =
    pmax(df$FW_Flag[df$row_id %in% df06$row_id], df06$FW_Flag)
  
  df06$Outlier = "Good"
  df06$Outlier[df06$HW_Flag == 2] = "HW_Area_Outlier"
  df06$Outlier[df06$FW_Flag == 2] = "FW_Area_Outlier"
  df06$Outlier[df06$HW_Flag == 3] = "Regression_Outlier"
  df06$Outlier[df06$FW_Flag == 3] = "Regression_Outlier"
  df06$Outlier[df06$HW_Flag == 2 & df06$FW_Flag == 2] = "FW_HW_Area_Outlier" 
  
  p3 = ggplot(df06, aes(x = FW_Area , y = HW_Area)) +
    geom_point(aes(color = Outlier),
               size = 0.6, alpha = 0.6) +
    
    stat_ellipse(aes(group = 1, linetype = "99%"),
                 type = "norm",
                 level = 0.99,
                 color = rangla_panjab[9],
                 linewidth = 0.8) +
    stat_ellipse(aes(group = 1, linetype = "95%"),
                 type = "norm",
                 level = 0.95,
                 color = "deeppink4",
                 linewidth = 0.8) +
    stat_ellipse(aes(group = 1, linetype = "90%"),
                 type = "norm",
                 level = 0.90,
                 color = rangla_panjab[5],
                 linewidth = 0.8) +
    scale_color_manual(
      name = "Outlier Category",
      values = outlier_colors,
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = "Confidence ellipse",
      values = c(
        "90%" = "solid",
        "95%" = "dashed",
        "99%" = "dashed"
      )
    )+
    guides(
      color = guide_legend(
        title = "Outlier Category",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 2)
      ),
      linetype = guide_legend(
        title = "Confidence ellipse",
        override.aes = list(linewidth = 1.2)
      )
    )+
    xlab(expression(log[10](Forewing~area~(mm^2)))) +
    ylab(expression(log[10](Hindwing~area~(mm^2)))) +
    coord_fixed(xlim = c(0,3), ylim = c(0,3))+
    
    theme_bw()+
    theme(
      legend.position = "none",   # ← ADD THIS
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.key.size = unit(0.35, "cm")
    )
  plot = (p1 + p2 + p3) +
    plot_layout(ncol = 3)
  print(plot)
  
 
  
}