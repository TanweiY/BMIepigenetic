################### 1. remake the heatmap by only looking at 1-4 years ago and more than 5 years ##
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/spearmann_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/df_multyearvor_list.RData")

df_bmi14<-list()
summary(df_multyearvor_list[[1]]$year_diff)
for (i in 1:20){
  df_bmi14[[i]]<-subset(df_multyearvor_list[[i]],
                        df_multyearvor_list[[i]]$year_diff>=1 & df_multyearvor_list[[i]]$year_diff <=4)
  

}

df<-df_bmi14[[1]] ## 846


scores <- c("PI_mdson", "PI_mc", "PI_hmn", "PI_do", "PI_mzbach")
bmi_measure <- c("BMI_dyna")

# Define the lists of dataframes and corresponding time periods
dataframe_lists <- list(
  df_bmi14 = df_bmi14

)

time_periods <- c("1-4 y ago")

# Initialize an empty list to store results
all_results_yr <- list()

# Iterate through the dataframe lists and time periods
for (i in seq_along(dataframe_lists)) {
  # Get the current list of dataframes
  current_list <- dataframe_lists[[i]]
  
  # Get the corresponding time period
  time_period <- time_periods[i]
  
  # Calculate Spearman correlations for the current list
  result <- calculate_spearman(current_list, scores, bmi_measure)
  
  # Add the time period column
  result$time_period <- time_period
  
  # Append the result to the results list
  all_results_yr[[i]] <- result
}

# Combine all results into a single dataframe
final_resultsyr <- do.call(rbind, all_results_yr)

# Print the combined results
print(final_results)
final_resultsyr<-final_resultsyr[order(final_resultsyr$score, final_resultsyr$time_period), ]
save(final_resultsyr, 
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_14bdiagnosis.RData')

## start plotting the updated heatmap ##
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_allsample.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_14bdiagnosis.RData")

final_resultsyr$BMI_measure<-final_resultsyr$time_period
#final_resultsyr$BMI_measure<-gsub(" ago", "", final_resultsyr$BMI_measure)
final_resultsyr<-final_resultsyr[, 1:4]
result<-subset(result, BMI_measure ==  'BMI'|BMI_measure ==  'BMI_5_14earlier')

## bind all results together
hm_all<-rbind(result, final_resultsyr)
# start ploting 
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library("reshape") 
library(dplyr)
library(tidyr)
library(cgwtools)

hm_all$BMI_measure<-factor(hm_all$BMI_measure, 
                           levels = rev(c('BMI', '1-4 y ago', 'BMI_5_14earlier')),
                           labels = rev(c('At diagnosis', '1-4y ago', '5-14y ago')))

hm_all$score<-factor(hm_all$score,
                     levels = c("PI_mdson",  "PI_do", "PI_hmn",  "PI_mc",  "PI_mzbach"),
                     labels = c("mBMI-135", "mBMI-398","mBMI-400", "mBMI-1109", "mBMI-3506"))

options(digits = 3)
summary(hm_all$avg_p_value)
quartiles <- quantile(hm_all$avg_p_value, probs = c(0.25, 0.5, 0.75))
quartiles[1]
hm_all$avg_p_value <- round(hm_all$avg_p_value, 4)


ggplot(hm_all, aes(score, BMI_measure)) +                           
  geom_tile(aes(fill = avg_spearman_coeff), colour = "white") +
  scale_fill_gradient(high = "#CE2029", low = "#FFAAAA")+
  geom_text(aes(label = round(avg_spearman_coeff, 2)), size = 4, color = "black") +  # Add text labels
  labs(title = '', x='Methylation-based BMI score', y='BMI measurement', fill = "Spearman's coefficient")+
  scale_x_discrete(position = "top")+
  theme( axis.text = element_text(size = 12), axis.title=element_text(size=12),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=13, face = 'bold'),
         legend.position = "none",
         #axis.title.y = element_blank(),
         #axis.title.x= element_blank(),
         plot.margin = margin(0, 0, 0, 0, "cm"))

## 4 * 6.5

######################### 2. forest plot, add 5-14 year ago ########
library(forestplot)
library(readxl)

## forest plot for self-reported BMI
selfdf <- read_excel("methy_BMI/coxforest_selfreport.xlsx")

colnames(selfdf)

ds_data<-selfdf[, c(3:5)]
ds_data$HR<-as.numeric(ds_data$HR)
ds_data$LL<-round(as.numeric(ds_data$LL), digits = 2)
ds_data$UL<-round(as.numeric(ds_data$UL), digits = 2)

colnames(ds_data)<-c("mean", "lower", "upper")
ds_data<-ds_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(selfdf))
m1fp_text<-subset(selfdf, select = c("Measurement", "HRCI"))

m1fp_table <- cbind(
  c("Measurement",  m1fp_text$Measurement),
  c("HR (95%CI)", m1fp_text$HRCI))

forestplot(m1fp_table,
               ds_data, 
               is.summary = c(TRUE, TRUE, TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,3), 
                              TRUE, TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,3),
                              TRUE, TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,3)),
               boxsize = 0.4,
               align=c("l","c","l","c","c", "c","c"),
               graph.pos = 2,
               line.margin = unit(0.5, "cm"),
               colgap = unit(6, "mm"),
               col=fpColors(box="#D73027",
                            line="#4575B4"),
               hrzl_lines = list("2" = gpar(lty = 1)),
               clip= c(0.1, 2),
               xlab = "Hazard ratio",
               zero=1,
               graphwidth=  unit(5, "cm"),
               new_page = TRUE,
               lwd.ci=1,  
               lwd.xaxis=1,  
               lwd.zero=0.5,
               vertices = TRUE,
               ci.vertices.height= 0.2,
               txt_gp=fpTxtGp(label=gpar(fontsize=12),
                              ticks=gpar(fontsize=12),
                              xlab=gpar(fontsize=12)))

p1
# A4, protriat

#### add BMI 5-14 years ago in the dose-response curves #####

####
library(rms)
library(ggplot2)
library(plyr)
library(ggpubr)
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
df<-blood_all[[1]]
summary(df)
## randomly selected one dataset to plot the curves
dd <- datadist(df)
options(datadist="dd")

AIC <- Inf  # Initialize AIC with a high value
nk <- 3     # Default nk value

for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all) ~ rcs(BMI, i),
             data = df, x = TRUE)
  tmp <- AIC(fit)  # Use AIC function from rms package
  if (tmp < AIC) {
    AIC <- tmp
    nk <- i
  }
}

print(nk)
dput(names(df))
### plot the association with CSC
df$TNM_adj <- as.factor(as.character(df$TNM_adj))
is.ordered(df$TNM_adj)
cont_var <- c("BMI", "BMI_5_14earlier", "PI_mdson",  "PI_mc", "PI_hmn", "PI_do",  "PI_mzbach")
adjust_covariates <- "+Age_diag + Sex + Schooling_years + TNM_adj + Location + chemradther +
  neotreat+ smoking+physical_activity_lifeav + lifetimealcohol_day + cardiovad + NSAIDS +
  statins + diabetesT2 + high_blood_pressure + hyperlipidemia+other_cancer"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each BMI measure
for (var in cont_var) {
  
  # Fit the model
  formula <- as.formula(paste("Surv(time = time_blood, time2 = timeD, death_all) ~ rcs(", var, ", 3)", adjust_covariates))
  fit <- cph(formula, data = df, x = TRUE)
  
  # Generate predictions
  surv2 <- Predict(fit, name = var, ref.zero = TRUE, fun = exp)
  
  # Calculate quartiles for the BMI variable
  t <- quantile(df[[var]], probs = seq(0, 1, 1/4))[c(2, 3, 4)]
  
  # Create the ggplot for the current BMI measure
  plot <- ggplot(surv2,adj.subtitle = FALSE, aes(x = !!sym(var), y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "pink") +
    geom_line(color = "red") +
    labs(x = var, y = "Hazard ratio", title = paste(var, "(OS)")) +
   # geom_vline(xintercept = t, linetype = "dashed", color = "gray") +
    coord_cartesian(ylim = c(0, 2.5)) +  # Use this instead of ylim()
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "top",
          legend.spacing.x = unit(1.2, 'cm'),
          legend.text = element_text(size = 11),
           axis.line.x = element_line(color = "black", size = 0.8),
      axis.line.y = element_line(color = "black", size = 0.8),
       axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
      axis.title.x = element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")))
    #+annotate("text", x = t[1], y = 0, label = "Q1", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    #annotate("text", x = t[2], y = 0, label = "Q2", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    #annotate("text", x = t[3], y = 0, label = "Q3", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic")
  
  # Add the plot to the list
  plot_list[[var]] <- plot
}

plot_list[[3]]

os<-ggarrange(plot_list[[1]], plot_list[[2]],plot_list[[3]],
              plot_list[[4]], plot_list[[5]],plot_list[[6]], plot_list[[7]], nrow = 7, ncol = 1)

## 6 *11 

##### cancer-specific survival ###
df$TNM_adj <- as.factor(as.character(df$TNM_adj))
is.ordered(df$TNM_adj)
cont_var <- c("BMI", "BMI_5_14earlier", "PI_mdson",  "PI_mc", "PI_hmn", "PI_do",  "PI_mzbach")
adjust_covariates <- "+Age_diag + Sex + Schooling_years + TNM_adj + Location + chemradther +
  neotreat+ smoking+physical_activity_lifeav + lifetimealcohol_day + cardiovad + NSAIDS +
  statins + diabetesT2 + high_blood_pressure + hyperlipidemia+other_cancer"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each BMI measure
for (var in cont_var) {
  
  # Fit the model
  formula <- as.formula(paste("Surv(time = time_blood, time2 = timeD, death_crccp==1) ~ rcs(", var, ", 3)", adjust_covariates))
  fit <- cph(formula, data = df, x = TRUE)
  
  # Generate predictions
  surv2 <- Predict(fit, name = var, ref.zero = TRUE, fun = exp)
  
  # Calculate quartiles for the BMI variable
  t <- quantile(df[[var]], probs = seq(0, 1, 1/4))[c(2, 3, 4)]
  
  # Create the ggplot for the current BMI measure
  plot <- ggplot(surv2,adj.subtitle = FALSE, aes(x = !!sym(var), y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "pink") +
    geom_line(color = "red") +
    labs(x = var, y = "Hazard ratio", title = paste(var, "(CSS)")) +
    #geom_vline(xintercept = t, linetype = "dashed", color = "gray") +
    coord_cartesian(ylim = c(0, 2.5)) +  # Use this instead of ylim()
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "top",
          legend.spacing.x = unit(1.2, 'cm'),
          legend.text = element_text(size = 11),
           axis.line.x = element_line(color = "black", size = 0.8),
            axis.line.y = element_line(color = "black", size = 0.8),
           axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm"))) 
   # +annotate("text", x = t[1], y = 0, label = "Q1", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
   # annotate("text", x = t[2], y = 0, label = "Q2", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
   # annotate("text", x = t[3], y = 0, label = "Q3", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic")
  
  # Add the plot to the list
  plot_list[[var]] <- plot
}


css<-ggarrange(plot_list[[1]], plot_list[[2]],plot_list[[3]],
               plot_list[[4]], plot_list[[5]],plot_list[[6]],plot_list[[7]],
               nrow = 7, ncol = 1)

drplot<-ggarrange(os, css, ncol  = 2)

drplot





























