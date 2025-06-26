load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_allv20.RData")
## only use one dataset, because the percentages of missing values are not too high
df<-blood_all[[1]]
summary(df)
dput(names(df))
# Define the scores and BMI measures
scores <- c("PI_mdson", "PI_mc", "PI_hmn", "PI_do", "PI_mzbach")
bmi_measure <- c("BMI", "BMI_5_14earlier", "WYO_td", "WYO_tr5yr", 
                 "medianBMI_tdiag", "meanBMI_tdiag", "medianBMI_t514", "meanBMI_t514")

# standardization
df[scores]<-scale(df[scores])
df[bmi_measure]<-scale(df[bmi_measure])
summary(df)


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

write_xlsx(adjusted_results,
           '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/results/linearegression_adjusted.xlsx')


