source("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)
library(survival)

load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_allv20.RData")
## only use one dataset, because the percentages of missing values are not too high
df<-blood_all[[9]]
summary(df)
dput(names(df))

df<-within.data.frame(df, {
  BMI_cat<- cut(BMI, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  BMI5_14er_cat<- cut(BMI_5_14earlier, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  meanBMI_tdiagcat<- cut(meanBMI_tdiag, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  meanBMI_t514cat<- cut(meanBMI_t514, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  
  WYO_tdqt<-factor(ntile(WYO_td, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  WYO_tr5yrqt<-factor(ntile(WYO_tr5yr, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  PI_mdsonqt<-factor(ntile(PI_mdson, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  PI_mcqt<-factor(ntile(PI_mc, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  PI_hmnqt<-factor(ntile(PI_hmn, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  PI_doqt<-factor(ntile(PI_do, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  PI_mzbachqt<-factor(ntile(PI_mzbach, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))

})

summary(df)
save(df,
     file = '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_progdf.RData')

# Define the scores and BMI measures
cont_var<-c("BMI", "BMI_5_14earlier", "WYO_td","WYO_tr5yr","meanBMI_tdiag","meanBMI_t514", "PI_mc", 
            "PI_hmn", "PI_do", "PI_mdson", "PI_mzbach")
  
df[cont_var]<-scale(df[cont_var])


scores <- c("BMI", "BMI_cat", "BMI_5_14earlier", "BMI5_14er_cat", "WYO_td", "WYO_tdqt",  "WYO_tr5yr", "WYO_tr5yrqt",
                 "meanBMI_tdiag", "meanBMI_tdiagcat", "meanBMI_t514", "meanBMI_t514cat",
                 "PI_mdson", "PI_mdsonqt", "PI_mc", "PI_mcqt", "PI_hmn","PI_hmnqt", "PI_do","PI_doqt",  "PI_mzbach", "PI_mzbachqt")

outcomes<-c('OS', 'DOC', 'CSS', 'TTR')

blood_cox<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxfull(data = df, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox<- rbind(blood_cox, result)
  }
}  

# Replace '[' with '(' and ']' with ')' in the dataframe
blood_cox <- as.data.frame(lapply(blood_cox, function(x) {
  if (is.character(x)) {
    x <- gsub("\\[", "(", x)  # Replace '[' with '('
    x <- gsub("\\]", ")", x)  # Replace ']' with ')'
    x
  } else {
    x  # Leave other data types unchanged
  }
}), stringsAsFactors = FALSE)

# Reshape from long to wide format
result_wide <- reshape(blood_cox, 
                   timevar = "outcome",    # Column that defines new wide columns
                   idvar = "score",        # Identifier column (stays as rows)
                   direction = "wide",     # Reshape to wide format
                   v.names = c("aHR.95..CI.", "aPvalue") )

write_xlsx(result_wide,
          '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/results/cox_results_all.xlsx')

setwd('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform')

######################## methylation-based scores ##################
cont_var <- c("PI_mc", "PI_hmn", "PI_do", "PI_mdson", "PI_mzbach")
adjust_covariates <- "+Age_diag + Sex + Schooling_years + TNM_adj + crc2sites + chemradther + smoking +
                      physical_activity_lifeav + hormone_replace + cardiovad + NSAIDS + statins + other_cancer"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each BMI measure
for (var in cont_var) {
  
  # Fit the model
  formula <- as.formula(paste("Surv(timeD, death_all) ~ rcs(", var, ", 3)", adjust_covariates))
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
    geom_vline(xintercept = t, linetype = "dashed", color = "gray") +
    ylim(0, 2.5) +
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
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm"))) +
    annotate("text", x = t[1], y = 0, label = "Q1", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    annotate("text", x = t[2], y = 0, label = "Q2", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    annotate("text", x = t[3], y = 0, label = "Q3", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic")
  
  # Add the plot to the list
  plot_list[[var]] <- plot
}

ggarrange(plot_list[[1]], plot_list[[2]],plot_list[[3]],
          plot_list[[4]], plot_list[[5]])












  