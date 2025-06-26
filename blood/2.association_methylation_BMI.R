load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
########## whole sampe ################
# Define parameters
scores <- c("PI_mdson", "PI_mc", "PI_hmn", "PI_do", "PI_mzbach")
bmi_measure <- c("BMI", "BMI_5_14earlier", "meanBMI_tdiag",  "meanBMI_t514")

# Initialize a list to store results from all datasets
all_results <- list()

# Loop through each imputed dataset
for (i in seq_along(blood_all)) {
  df <- blood_all[[i]]
  
  # Standardize scores and BMI measures
  df[scores] <- scale(df[scores])
  df[bmi_measure] <- scale(df[bmi_measure])
  
  # Calculate Spearman correlations for this dataset
  results <- list()
  for (score in scores) {
    for (bmi in bmi_measure) {
      corr_test <- cor.test(df[[score]], df[[bmi]], method = "spearman")
      results <- append(results, list(data.frame(
        score = score,
        BMI_measure = bmi,
        spearman_coeff = corr_test$estimate,
        p_value = corr_test$p.value
      )))
    }
  }
  
  # Combine results for this dataset
  dataset_results <- do.call(rbind, results)
  dataset_results$dataset <- i  # Add dataset identifier
  all_results[[i]] <- dataset_results
}

# Combine results across all datasets
combined_results <- do.call(rbind, all_results)

# Average the coefficients and p-values across datasets
averaged_results <- combined_results %>%
  group_by(score, BMI_measure) %>%
  summarise(
    avg_spearman_coeff = mean(spearman_coeff, na.rm = TRUE),
    avg_p_value = mean(p_value, na.rm = TRUE),
    .groups = "drop"
  )

averaged_results$avg_p_value<-NULL # p all <0.0005
averaged_results$avg_spearman_coeff<-round(averaged_results$avg_spearman_coeff, 2)








# Save the result to an Excel file
write_xlsx(averaged_results, "averaged_spearman_results.xlsx")
########


# Print the result
print(wide_df)
wide_df <- wide_df %>%
  mutate(score = factor(score, levels = scores)) %>%
  arrange(score)

wide_df$Study<-c('Mendelson et al. 2017', 'McCartney et al. 2018', 'Hamilton et al. 2019', 
                 'Do et al. 2022', 'Merzbacher et al. 2023')
wide_df$Markers<-c('mBMI-135','mBMI-1109', 'mBMI-400', 'mBMI-398', 'mBMI-3506')

save(wide_df,
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_allsample.RData')

################ wanna make a heatmap out of this #####




############ linear regression #############

# Initialize an empty list to store the results
results <- list()

# Nested loop to run regression for each score and BMI measure
for (score in scores) {
  for (bmi in bmi_measure) {
    # Fit the linear regression model
    formula <- as.formula(paste(bmi, "~", score, 
                                "+ Age_diag + Sex + Schooling_years + TNM_adj + crc2sites + chemradther + smoking+physical_activity_lifeav+ hormone_replace+cardiovad+NSAIDS+statins+other_cancer"))
    model <- lm(formula, data = df)
    
    # Extract regression results
    estimate <- summary(model)$coefficients[score, "Estimate"]
    std_error <- summary(model)$coefficients[score, "Std. Error"]
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    p_value <- summary(model)$coefficients[score, "Pr(>|t|)"]
    r_squared <- summary(model)$r.squared
    
    # Create a string for the estimate with 95% CI
    estimate_ci <- paste0(round(estimate, 2), " (", 
                          round(ci_lower, 2), ", ", 
                          round(ci_upper, 2), ")")
    
    # Append results as a data frame to the list
    results <- append(results, list(data.frame(
      score = score,
      BMI_measure = bmi,
      estimate_CI = estimate_ci,
      p_value = round(p_value, 3),
      R2 = round(r_squared, 3)
    )))
  }
}

# Combine all results into a single data frame
adjusted_results <- do.call(rbind, results)
adjusted_results <- adjusted_results %>%
  pivot_wider(
    names_from = BMI_measure,  # Columns to create new column names
    values_from = c(estimate_CI, p_value, R2)  # Values to fill in the cells
  )

adjusted_results<-adjusted_results[, c(1:9)]
adjusted_results$measure<-'linear regression'
colnames(adjusted_results)<-colnames(wide_df)
wide_df[, c(2:9)]<-lapply(wide_df[, c(2:9)], as.character)

myselfasso_allsample<-bind_rows(wide_df, adjusted_results)

write_xlsx(myselfasso_allsample, 
           '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/results/myselfasso_allsample.xlsx')


########################### stratified analyses by data collection 
load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_allv20.RData")

########## whole sampe ################
df<-blood_all[[1]]
summary(df$time_blood)
df$time_blood<-df$time_blood/30

scores<-c("PI_mdson", "PI_mc", "PI_hmn", "PI_do","PI_mzbach")

bmi_measure<-c("BMI", "BMI_5_14earlier", 
               "WYO_td", "WYO_tr5yr", "medianBMI_tdiag", "meanBMI_tdiag", "medianBMI_t514", 
               "meanBMI_t514")

# standardization
df[scores]<-scale(df[scores])
df[bmi_measure]<-scale(df[bmi_measure])


df_e<-subset(df, time_blood<3) # 1233

# Nested loop to calculate Spearman correlations
results <- list()

for (score in scores) {
  for (bmi in bmi_measure) {
    # Calculate Spearman correlation
    corr_test <- cor.test(df_e[[score]], df_e[[bmi]], method = "spearman")
    
    # Extract results and append to list
    results <- append(results, list(data.frame(
      score = score,
      BMI_measure = bmi,
      spearman_coeff = corr_test$estimate,
      p_value = corr_test$p.value
    )))
  }
}

# Combine all results into a single data frame
coeff_all <- do.call(rbind, results)
coeff_all$spearman_coeff<-round(coeff_all$spearman_coeff, 2)
coeff_all$p_value<-round(coeff_all$p_value, 3)
coeff_all$p_value[coeff_all$p_value ==0]<-'<0.0001'
coeff_all$p_value<-NULL
library(tidyr)
library(writexl)
wide_df <- coeff_all %>%
  pivot_wider(
    names_from = BMI_measure,  # Columns to create new column names
    values_from = c(spearman_coeff)  # Values to fill in the cells
  )

wide_df$measure<-'spearmann'
wide_df$group<-'Blood<3months'

#### linear regression #####
# Initialize an empty list to store the results
results <- list()

# Nested loop to run regression for each score and BMI measure
for (score in scores) {
  for (bmi in bmi_measure) {
    # Fit the linear regression model
    formula <- as.formula(paste(bmi, "~", score, 
                                "+ Age_diag + Sex + Schooling_years + TNM_adj + crc2sites + chemradther + smoking+physical_activity_lifeav+ hormone_replace+cardiovad+NSAIDS+statins+other_cancer"))
    model <- lm(formula, data = df_e)
    
    # Extract regression results
    estimate <- summary(model)$coefficients[score, "Estimate"]
    std_error <- summary(model)$coefficients[score, "Std. Error"]
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    p_value <- summary(model)$coefficients[score, "Pr(>|t|)"]
    r_squared <- summary(model)$r.squared
    
    # Create a string for the estimate with 95% CI
    estimate_ci <- paste0(round(estimate, 2), " (", 
                          round(ci_lower, 2), ", ", 
                          round(ci_upper, 2), ")")
    
    # Append results as a data frame to the list
    results <- append(results, list(data.frame(
      score = score,
      BMI_measure = bmi,
      estimate_CI = estimate_ci,
      p_value = round(p_value, 3),
      R2 = round(r_squared, 3)
    )))
  }
}

# Combine all results into a single data frame
adjusted_results <- do.call(rbind, results)
adjusted_results <- adjusted_results %>%
  pivot_wider(
    names_from = BMI_measure,  # Columns to create new column names
    values_from = c(estimate_CI, p_value, R2)  # Values to fill in the cells
  )

adjusted_results<-adjusted_results[, c(1:9)]
adjusted_results$measure<-'linear regression'
adjusted_results$group<-'Blood<3months'

colnames(adjusted_results)<-colnames(wide_df)
wide_df[, c(2:9)]<-lapply(wide_df[, c(2:9)], as.character)

myselfasso_BCL3M<-bind_rows(wide_df, adjusted_results)
write_xlsx(myselfasso_BCL3M, 
           '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/results/myselfasso_BCL3M.xlsx')



load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_allv20.RData")

df<-blood_all[[1]]
summary(df$time_blood)
df$time_blood<-df$time_blood/30

scores<-c("PI_mdson", "PI_mc", "PI_hmn", "PI_do","PI_mzbach")

bmi_measure<-c("BMI", "BMI_5_14earlier", 
               "WYO_td", "WYO_tr5yr", "medianBMI_tdiag", "meanBMI_tdiag", "medianBMI_t514", 
               "meanBMI_t514")

# standardization
df[scores]<-scale(df[scores])
df[bmi_measure]<-scale(df[bmi_measure])


df_e<-subset(df, time_blood>=3) # 893

# Nested loop to calculate Spearman correlations
results <- list()

for (score in scores) {
  for (bmi in bmi_measure) {
    # Calculate Spearman correlation
    corr_test <- cor.test(df_e[[score]], df_e[[bmi]], method = "spearman")
    
    # Extract results and append to list
    results <- append(results, list(data.frame(
      score = score,
      BMI_measure = bmi,
      spearman_coeff = corr_test$estimate,
      p_value = corr_test$p.value
    )))
  }
}

# Combine all results into a single data frame
coeff_all <- do.call(rbind, results)
coeff_all$spearman_coeff<-round(coeff_all$spearman_coeff, 2)
coeff_all$p_value<-round(coeff_all$p_value, 3)
coeff_all$p_value[coeff_all$p_value ==0]<-'<0.0001'
coeff_all$p_value<-NULL
library(tidyr)
library(writexl)
wide_df <- coeff_all %>%
  pivot_wider(
    names_from = BMI_measure,  # Columns to create new column names
    values_from = c(spearman_coeff)  # Values to fill in the cells
  )

wide_df$measure<-'spearmann'
wide_df$group<-'Blood>3months'

#### linear regression #####
# Initialize an empty list to store the results
results <- list()

# Nested loop to run regression for each score and BMI measure
for (score in scores) {
  for (bmi in bmi_measure) {
    # Fit the linear regression model
    formula <- as.formula(paste(bmi, "~", score, 
                                "+ Age_diag + Sex + Schooling_years + TNM_adj + crc2sites + chemradther + smoking+physical_activity_lifeav+ hormone_replace+cardiovad+NSAIDS+statins+other_cancer"))
    model <- lm(formula, data = df_e)
    
    # Extract regression results
    estimate <- summary(model)$coefficients[score, "Estimate"]
    std_error <- summary(model)$coefficients[score, "Std. Error"]
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    p_value <- summary(model)$coefficients[score, "Pr(>|t|)"]
    r_squared <- summary(model)$r.squared
    
    # Create a string for the estimate with 95% CI
    estimate_ci <- paste0(round(estimate, 2), " (", 
                          round(ci_lower, 2), ", ", 
                          round(ci_upper, 2), ")")
    
    # Append results as a data frame to the list
    results <- append(results, list(data.frame(
      score = score,
      BMI_measure = bmi,
      estimate_CI = estimate_ci,
      p_value = round(p_value, 3),
      R2 = round(r_squared, 3)
    )))
  }
}

# Combine all results into a single data frame
adjusted_results <- do.call(rbind, results)
adjusted_results <- adjusted_results %>%
  pivot_wider(
    names_from = BMI_measure,  # Columns to create new column names
    values_from = c(estimate_CI, p_value, R2)  # Values to fill in the cells
  )

adjusted_results<-adjusted_results[, c(1:9)]
adjusted_results$measure<-'linear regression'
adjusted_results$group<-'Blood>3months'

colnames(adjusted_results)<-colnames(wide_df)
wide_df[, c(2:9)]<-lapply(wide_df[, c(2:9)], as.character)

myselfasso_BCM3M<-bind_rows(wide_df, adjusted_results)
write_xlsx(myselfasso_BCM3M, 
           '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/results/myselfasso_BCM3M.xlsx')


















