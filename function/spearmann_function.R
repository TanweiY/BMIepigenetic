library(dplyr)

# Define the function
calculate_spearman <- function(data_list, scores, bmi_measure) {
  # Initialize a list to store results from all datasets
  all_results <- list()
  
  # Loop through each dataset in the list
  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    
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
  
  # Optional: Drop p_value column and round coefficients
  #averaged_results$avg_p_value <- NULL  # p-values are assumed < 0.0005
  averaged_results$avg_spearman_coeff <- round(averaged_results$avg_spearman_coeff, 2)
  
  return(averaged_results)
  
  
}
