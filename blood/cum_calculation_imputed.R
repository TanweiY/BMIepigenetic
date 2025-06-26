load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_BMIcum20raw.RData")

# Initialize an empty list to store results
cumBMI_result_list <- list()

# Loop through each dataset in BMIcum
for (dataset_index in 1:20) {
  # Load the current dataset (e.g., BMIcum[[1]], BMIcum[[2]], ..., BMIcum[[20]])
  df <- BMIcum[[dataset_index]]
  
  # Convert from wide to long format
  long_data <- reshape(
    df,
    varying = list(c("BMI_20", "BMI_30", "BMI_40", "BMI_50", "BMI_60", "BMI_70", "BMI_80")), 
    v.names = "BMI_dyna",
    timevar = "Age_measure",
    times = c(20, 30, 40, 50, 60, 70, 80),
    direction = "long"
  )
  
  long_data <- long_data[order(long_data$tn_id), ]
  long_data <- long_data %>% 
    group_by(tn_id) %>%
    filter(Age_measure <= Age_diag)
  
  # Initialize an empty list to store results for each dataset
  result_list <- list()
  
  # Get unique 'id' values
  unique_ids <- unique(long_data$id)
  
  # Iterate over each 'id'
  for (i in unique_ids) {
    subset_df <- subset(long_data, id == i)
    
    # Perform linear interpolation
    interpolation <- approx(
      x = subset_df$Age_measure,
      y = subset_df$BMI_dyna,
      xout = seq(20, unique(subset_df$Age_diag)),
      rule = 2
    )
    
    group_result <- data.frame(
      tn_id = unique(subset_df$tn_id),
      Age_diag = unique(subset_df$Age_diag),
      Year_recent = unique(subset_df$Year_recent),
      Age_measure = interpolation$x,
      BMI_dyna = interpolation$y
    )
    
    # Append the result to the list for this dataset
    result_list[[i]] <- group_result
  }
  
  # Combine all results for this dataset into a single dataframe
  final_result <- do.call(rbind, result_list)
  final_result$eBMI_dyn <- final_result$BMI_dyna - 25
  final_result$eBMI_dyn[final_result$eBMI_dyn < 0] <- 0
  
  # Calculate WYO, cuBMI, etc.
  
  # 1. WYO_tdiagnosis: eBMI from 20 to diagnosis
  WYO_td <- aggregate(eBMI_dyn ~ tn_id, 
                      data = final_result, FUN = sum, na.rm = TRUE)
  colnames(WYO_td)[2] <- 'WYO_td'
  
  # 2. WYO_tr5yr: eBMI from 20 to 5-14 years before diagnosis
  final_result$lower_tens <- floor(final_result$Year_recent / 10) * 10
  BMI_514 <- final_result %>% 
    group_by(tn_id) %>%
    filter(Age_measure <= lower_tens)
  
  WYO_tr5yr <- aggregate(eBMI_dyn ~ tn_id, 
                         data = BMI_514, FUN = sum, na.rm = TRUE)
  colnames(WYO_tr5yr)[2] <- 'WYO_tr5yr'
  
  # 3. Median and Mean BMI from 20 to diagnosis / 5-14 years before diagnosis
  cuBMI_mediandiagn <- aggregate(BMI_dyna ~ tn_id, 
                                 data = final_result, FUN = median, na.rm = TRUE)
  
  cuBMI_meandiagn <- aggregate(BMI_dyna ~ tn_id, 
                               data = final_result, FUN = mean, na.rm = TRUE)
  
  cuBMI_medianearl <- aggregate(BMI_dyna ~ tn_id, 
                                data = BMI_514, FUN = median, na.rm = TRUE)
  
  cuBMI_meanearl <- aggregate(BMI_dyna ~ tn_id, 
                              data = BMI_514, FUN = mean, na.rm = TRUE)
  
  colnames(cuBMI_meandiagn)[2] <- 'meanBMI_tdiag'
  colnames(cuBMI_mediandiagn)[2] <- 'medianBMI_tdiag'
  colnames(cuBMI_medianearl)[2] <- 'medianBMI_t514'
  colnames(cuBMI_meanearl)[2] <- 'meanBMI_t514'
  
  # Bind all results together for this dataset
  cumBMI_result <- merge(WYO_td, WYO_tr5yr, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_mediandiagn, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_meandiagn, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_medianearl, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_meanearl, by = 'tn_id')
  
  # Store the result of this dataset in the list
  cumBMI_result_list[[dataset_index]] <- cumBMI_result
}

d<-cumBMI_result_list[[1]]
summary(d)
save(cumBMI_result_list, 
     file = "/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/cumBMI_result_list20.RData")








