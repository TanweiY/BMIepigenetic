library(survival)

# Modified multicox function to extract HR and SE only for the 'score' variable (categorical or continuous)
full_multicox <- function(data, entrytime, score, outcome){
  if (!outcome %in% c('OS', 'DOC', 'CSS')){
    stop("Invalid outcome. Please choose 'OS', 'DOC', or 'CSS'.")
  }
  
  # Build formula (with confounders, but we'll only extract 'score')
  fml_base <- paste0("Surv(time = ", entrytime, ", time2 = timeD, ",
                    if (outcome == 'OS') "death_all)" else
                      if (outcome == 'DOC') "death_crccp==2)" else "death_crccp==1)", 
                    "~Age_diag + Sex + Schooling_years + TNM_adj + Location + chemradther +
                    neotreat+ smoking+physical_activity_lifeav + lifetimealcohol_day + cardiovad + NSAIDS +
                    statins + diabetesT2 + high_blood_pressure + hyperlipidemia+other_cancer+ ")
  
  fml <- as.formula(paste0(fml_base, score))
  
  # Fit Cox model
  model <- coxph(fml, data = data)
  
  # Extract HR and SE for only the 'score' variable
  model_summary <- summary(model)
  coef_table <- model_summary$coefficients
  
  # Filter only the rows corresponding to the 'score' variable
  # This will capture the HR and SE for all levels of a categorical variable or the continuous score
  score_rows <- grep(score, rownames(coef_table))
  
  # Extract HR and SE for the 'score'
  hr_se_table <- data.frame(
    HR = exp(coef_table[score_rows, "coef"]),
    SE = coef_table[score_rows, "se(coef)"],
    levels = rownames(coef_table)[score_rows]
  )
  
  return(hr_se_table)
}

# Function to run Cox models across dataframes and aggregate the results for the 'score' variable only
run_fullmulticox_impu <- function(dfs, entrytime, score, outcomes){
  results_list <- list()  # To store results for all score/outcome combinations
  
  for (outcome in outcomes) {
    all_results <- NULL  # Reset for each outcome
    
    for (df in dfs) {
      
      # Run the multicox function for each dataset
      result <- full_multicox(data = df, entrytime = entrytime, score = score, outcome = outcome)
      
      # Store HR and SE for all levels of the categorical/continuous 'score' variable across all datasets
      if (is.null(all_results)) {
        all_results <- result
      } else {
        all_results <- rbind(all_results, result)
      }
    }
    
    # Calculate mean HR and SE for each level of the 'score' across the datasets
    aggregated_results <- aggregate(all_results[, c("HR", "SE")], by = list(Level = all_results$levels), FUN = mean)
    
    # Calculate 95% CI for each level
    aggregated_results$Lower_95_CI <- exp(log(aggregated_results$HR) - 1.96 * aggregated_results$SE)
    aggregated_results$Upper_95_CI <- exp(log(aggregated_results$HR) + 1.96 * aggregated_results$SE)
    
    # Format HR (95% CI) with 2 decimal places
    aggregated_results$HR_95CI <- sprintf("%.2f (%.2f, %.2f)", aggregated_results$HR, aggregated_results$Lower_95_CI, aggregated_results$Upper_95_CI)
    
    # Add outcome column for clarity
    aggregated_results$Outcome <- outcome
    
    # Store the result
    results_list[[outcome]] <- aggregated_results
  }
  
  # Combine results into a single dataframe
  final_results <- do.call(rbind, results_list)
  return(final_results)
}

