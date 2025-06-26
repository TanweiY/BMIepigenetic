library(cgwtools)
library(dplyr)
library(writexl)
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/Cox_imputed_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/df_multyearvor_list.RData")

for (i in 1:20) {
  df_multyearvor_list[[i]]$BMI_dynacat<-cut(df_multyearvor_list[[i]]$BMI_dyna, 
                                breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  
  df_multyearvor_list[[i]]$BMI_dynacat<-relevel(df_multyearvor_list[[i]]$BMI_dynacat,
                                  ref = '20-25')
  
  
  df_multyearvor_list[[i]]$BMI_dynp5<-df_multyearvor_list[[i]]$BMI_dyna/5

}

save(df_multyearvor_list, 
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/df_multyearvor_list.RData')

df_bmi02<-list()
df_bmi24<-list()
df_bmi46<-list()
df_bmi68<-list()
df_bmi810<-list()
df_bmi1012<-list()

for (i in 1:20){
  df_bmi02[[i]]<-subset(df_multyearvor_list[[i]], df_multyearvor_list[[i]]$year_diff>0 & df_multyearvor_list[[i]]$year_diff <=2)
  df_bmi24[[i]]<-subset(df_multyearvor_list[[i]], df_multyearvor_list[[i]]$year_diff>2 & df_multyearvor_list[[i]]$year_diff <=4)
  df_bmi46[[i]]<-subset(df_multyearvor_list[[i]], df_multyearvor_list[[i]]$year_diff>4 & df_multyearvor_list[[i]]$year_diff <=6)
  df_bmi68[[i]]<-subset(df_multyearvor_list[[i]], df_multyearvor_list[[i]]$year_diff>6 & df_multyearvor_list[[i]]$year_diff <=8)
  df_bmi810[[i]]<-subset(df_multyearvor_list[[i]], df_multyearvor_list[[i]]$year_diff>8 & df_multyearvor_list[[i]]$year_diff <=10)
  df_bmi1012[[i]]<-subset(df_multyearvor_list[[i]], df_multyearvor_list[[i]]$year_diff>10 & df_multyearvor_list[[i]]$year_diff <=12)
}

# Define the scores and outcomes
scores <- c("BMI_dynp5",  "BMI_dynacat")
outcomes <- c('OS', 'DOC', 'CSS')

# Define the labels for time periods
time_periods <- c("0-2 y ago", "2-4 y ago", "4-6 y ago", "6-8 y ago", "8-10 y ago", "10-12 y ago")

# Initialize the final results data frame
all_multicox_results <- data.frame()

# Iterate through the datasets and time periods
for (time_idx in seq_along(time_periods)) {
  # Select the corresponding dataset for the current time period
  dfs <- switch(time_idx,
                df_bmi02,
                df_bmi24,
                df_bmi46,
                df_bmi68,
                df_bmi810,
                df_bmi1012)
  
  # Initialize a data frame for storing results of this time period
  multicox_results <- data.frame()
  
  # Run the multivariate Cox analysis for each score and outcome
  for (score in scores) {
    for (outcome in outcomes) {
      result <- run_fullmulticox_impu(dfs,  entrytime = 'time_blood', score = score, outcome = outcome)
      multicox_results <- rbind(multicox_results, result)
    }
  }
  
  # Process the results: order and reshape
  multicox_results <- multicox_results[order(multicox_results$Outcome), ]
  multicox_results <- multicox_results[, c(1, 6, 7)]
  multicox_results <- reshape(
    multicox_results,
    timevar = "Outcome",
    idvar = "Level",
    direction = "wide",
    v.names = c("HR_95CI"),
  )
  
  # Add the time period marker
  multicox_results$Time_Period <- time_periods[time_idx]
  multicox_results$Sample_size<-nrow(dfs[[1]])
  
  # Combine with the overall results
  all_multicox_results <- rbind(all_multicox_results, multicox_results)
}

# Final combined results
write_xlsx(all_multicox_results, 
           "C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/multicox_diagper2y.xlsx")



