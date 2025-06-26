source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/spearmann_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
scores <- c("PI_mdson", "PI_mc", "PI_hmn", "PI_do", "PI_mzbach")
bmi_measure <- c("BMI", "BMI_5_14earlier", "meanBMI_tdiag",  "meanBMI_t514")
result <- calculate_spearman(blood_all, scores, bmi_measure)

save(result, 
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_allsample.RData')

### calculate the association BMI by years before diagnosis
### first preprocess the data
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_BMIcum20raw.RData")

# Loop through each dataset in BMIcum
for (i in 1:20) {
  # Load the current dataset (e.g., BMIcum[[1]], BMIcum[[2]], ..., BMIcum[[20]])
  # Convert from wide to long format
  BMIcum[[i]] <- reshape(
    BMIcum[[i]],
    varying = list(c("BMI_20", "BMI_30", "BMI_40", "BMI_50", "BMI_60", "BMI_70", "BMI_80")), 
    v.names = "BMI_dyna",
    timevar = "Age_measure",
    times = c(20, 30, 40, 50, 60, 70, 80),
    direction = "long"
  )
  
  BMIcum[[i]] <- BMIcum[[i]][order(BMIcum[[i]]$tn_id), ]
  BMIcum[[i]] <- BMIcum[[i]] %>% 
    group_by(tn_id) %>%
    filter(Age_measure <= Age_diag)
  
  BMIcum[[i]]$year_diff<-BMIcum[[i]]$Age_diag-BMIcum[[i]]$Age_measure
  
}
  
save(BMIcum,
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_BMIcum20raw.RData")

### merge with clinical information 
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_imputatedpro.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_clin.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/scores.rdata")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_BMIcum20raw.RData")
df_clin<-subset(blood_clin, 
                select = c("tn_id", "lab_nummer", 
                           "late_days",  "time_blood",
                           "timeD", "timeM", "death_all", "death_crc", "timeD_recurr", "timeM_recurr", 
                           "recurr", "recurr_type", "recurr_cp", "death_crccp"))

colnames(ilm)[1]<-'lab_nummer'

df_clin<-merge(df_clin, ilm, by = 'lab_nummer', all.x = T)

df_multyearvor_list<-vector(20, mode = "list")

dput(names(BMIcum[[1]]))
for (i in 1:20) {
  # merge dateset
  BMIcum[[i]]<-subset(BMIcum[[i]], select = c("tn_id", "Age_diag", "Age_measure", "year_diff", "BMI_dyna"))

  blood_imputatedpro[[i]]<-subset(blood_imputatedpro[[i]],
                                  select =  c("tn_id", "Sex", "Schooling_years", "TNM_adj", "chemradther", 
                                              "neotreat", "Location", "lifetimealcohol_day", "smoking", "physical_activity_lifeav", 
                                              "hyperlipidemia", "statins", "NSAIDS", "diabetesT2", "high_blood_pressure", 
                                              "cardiovad", "other_cancer", "BMI", "BMI_5_14earlier", "height"))
  
  df_multyearvor_list[[i]]<-merge(blood_imputatedpro[[i]], df_clin, by = 'tn_id')
  df_multyearvor_list[[i]]<-merge(df_multyearvor_list[[i]], BMIcum[[i]], by = 'tn_id')
}

save(df_multyearvor_list,
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/df_multyearvor_list.RData")

## subset and iterate the function 
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/spearmann_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/df_multyearvor_list.RData")

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

scores <- c("PI_mdson", "PI_mc", "PI_hmn", "PI_do", "PI_mzbach")
bmi_measure <- c("BMI_dyna")

# Define the lists of dataframes and corresponding time periods
dataframe_lists <- list(
  df_bmi02 = df_bmi02,
  df_bmi24 = df_bmi24,
  df_bmi46 = df_bmi46,
  df_bmi68 = df_bmi68,
  df_bmi810 = df_bmi810,
  df_bmi1012 = df_bmi1012
)

time_periods <- c("0-2 y ago", "2-4 y ago", "4-6 y ago", "6-8 y ago", "8-10 y ago", "10-12 y ago")

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
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_year_diagnosis.RData')

## start plot the heat map
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_allsample.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/temp_results/spearmann_year_diagnosis.RData")
final_resultsyr$BMI_measure<-final_resultsyr$time_period
final_resultsyr$BMI_measure<-gsub(" ago", "", final_resultsyr$BMI_measure)
final_resultsyr<-final_resultsyr[, 1:4]
result<-subset(result, BMI_measure == 'BMI')

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
                           levels = rev(c('BMI', '0-2 y', '2-4 y', '4-6 y', '6-8 y', '8-10 y', '10-12 y')),
                           labels = rev(c('At diagnosis', '0-2y ago', '2-4y ago', '4-6y ago', '6-8y ago', '8-10y ago', '10-12y ago')))

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










  
  





