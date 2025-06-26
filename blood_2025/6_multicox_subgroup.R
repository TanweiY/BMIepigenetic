### subgroup analysis by blood collection time, stage (I/III vs IV; I/II vs. III/IV), gender, age, location (3 groups)
library(cgwtools)
library(dplyr)
library(writexl)
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/Cox_imputed_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")

for (i in 1:20){
  blood_all[[i]]$time_bloodm<-blood_all[[i]]$time_blood/30
}

bloodtl3<-list()
bloodtm3<-list()
stage12<-list()
stage34<-list()
stage13<-list()
stage4<-list()
male<-list()
female<-list()
agel70<-list()
agem70<-list()
discolon<-list()
proxcolon<-list()
rectum<-list()
df_bmil20<-list()
df_bmi2025<-list()
df_bmi2530<-list()
df_bmim30<-list()


for (i in 1:20){
  bloodtl3[[i]]<-subset(blood_all[[i]],blood_all[[i]]$time_bloodm <=3)
  bloodtm3[[i]]<-subset(blood_all[[i]],blood_all[[i]]$time_bloodm >3)
  stage12[[i]]<-subset(blood_all[[i]], blood_all[[i]]$TNM_adj =='I'|blood_all[[i]]$TNM_adj =='II')
  stage34[[i]]<-subset(blood_all[[i]], blood_all[[i]]$TNM_adj =='III'|blood_all[[i]]$TNM_adj =='IV')
  stage13[[i]]<-subset(blood_all[[i]], blood_all[[i]]$TNM_adj =='I'|blood_all[[i]]$TNM_adj =='II'|blood_all[[i]]$TNM_adj =='III')
  stage4[[i]]<-subset(blood_all[[i]], blood_all[[i]]$TNM_adj =='IV')
  male[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Sex =='Male')
  female[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Sex =='Female')
  agel70[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Age_diag <=70)
  agem70[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Age_diag >70)
  discolon[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Location == 'Distal colon')
  proxcolon[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Location == 'Proximal colon')
  rectum[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Location == 'Rectum')
  df_bmil20[[i]]<-subset(blood_all[[i]],blood_all[[i]]$BMI_cat =='<20')
  df_bmi2025[[i]]<-subset(blood_all[[i]],blood_all[[i]]$BMI_cat =='20-25')
  df_bmi2530[[i]]<-subset(blood_all[[i]],blood_all[[i]]$BMI_cat == '25-30')
  df_bmim30[[i]]<-subset(blood_all[[i]], blood_all[[i]]$BMI_cat == '≥30')
}

# Create a list to store all subgroup datasets
# Create a list of all subgroup lists
subgroup_lists <- list(
  bloodtl3 = bloodtl3,
  bloodtm3 = bloodtm3,
  stage12 = stage12,
  stage34 = stage34,
  stage13 = stage13,
  stage4 = stage4,
  male = male,
  female = female,
  agel70 = agel70,
  agem70 = agem70,
  #discolon = discolon, always have the problem: Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  :  exp overflow due to covariates
  proxcolon = proxcolon,
  rectum = rectum,
  df_bmil20 = df_bmil20,
  df_bmi2025 = df_bmi2025,
  df_bmi2530 = df_bmi2530,
  df_bmim30 = df_bmim30
)

## have to standardize all the scores once again
# Variables to standardize
variables_to_scale <- c("PI_mdson", "PI_mc", "PI_hmn", "PI_do", "PI_mzbach")

# Loop through each subgroup list
for (list_name in names(subgroup_lists)) {
  # Get the current subgroup list
  subgroup <- subgroup_lists[[list_name]]
  
  # Loop through each dataset in the subgroup list
  for (i in seq_along(subgroup)) {
    # Check if the dataset is not empty
    if (!is.null(subgroup[[i]]) && nrow(subgroup[[i]]) > 0) {
      # Loop through each variable to standardize
      for (var in variables_to_scale) {
        # Create the new variable name with "_sc" suffix
        new_var_name <- paste0(var, "_sc")
        # Standardize the variable and assign it to the new column
        subgroup[[i]][[new_var_name]] <- scale(subgroup[[i]][[var]])
      }
    }
  }
  
  # Update the original subgroup list with the standardized datasets
  subgroup_lists[[list_name]] <- subgroup
}

scores<-c("BMI_p5", "PI_mdson_sc",  "PI_mc_sc",  
          "PI_hmn_sc", "PI_do_sc",  "PI_mzbach_sc")

outcomes <- c('OS', 'DOC', 'CSS')

# Initialize an empty dataframe to store results
subgroup_multicox_results <- data.frame()

# Loop through each subgroup
for (subgroup_name in names(subgroup_lists)) {
  subgroup_data <- subgroup_lists[[subgroup_name]]
  
  # Loop through each score
  for (score in scores) {
    
    # Loop through each outcome
    for (outcome in outcomes) {
      
      # Run the Cox regression function for the current subgroup
      result <- run_fullmulticox_impu(
        dfs = subgroup_data,
        entrytime = 'time_blood',
        score = score,
        outcomes = outcome
      )
      
      # Add the subgroup information to the result
      result$Subgroup <- subgroup_name
      result$Score <- score  # Include score for clarity
      
      # Bind the result to the final dataframe
      subgroup_multicox_results <- rbind(subgroup_multicox_results, result)
    }
  }
}

# Reorder the results dataframe for better readability
subgroup_multicox_results <- subgroup_multicox_results[order(subgroup_multicox_results$Outcome), ]
subgroup_multicox_results<-subgroup_multicox_results[, c(1, 6, 7,8)]
multicox_results_sub<-reshape(subgroup_multicox_results, 
                          timevar = "Outcome",
                          idvar = c("Level", "Subgroup"),
                          direction = "wide",
                          v.names = c("HR_95CI"))

multicox_results_sub<-multicox_results_sub[order(multicox_results_sub$Level), ]

write_xlsx(multicox_results_sub, 
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/multicox_subgroup.xlsx')


### now deal with the sticking point ##
## randomly select one dataframe
library(cgwtools)
library(dplyr)
library(writexl)
library(tableone)
source("C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/Cox_complete_function.R")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
discolon<-list()
for (i in 1:20){
  discolon[[i]]<-subset(blood_all[[i]], blood_all[[i]]$Location == 'Distal colon')
  discolon[[i]]$PI_mdson_sc<-scale(discolon[[i]]$PI_mdson)
  discolon[[i]]$PI_mc_sc<-scale(discolon[[i]]$PI_mc)
  discolon[[i]]$PI_hmn_sc<-scale(discolon[[i]]$PI_hmn)
  discolon[[i]]$PI_do_sc<-scale(discolon[[i]]$PI_do)
  discolon[[i]]$PI_mzbach_sc<-scale(discolon[[i]]$PI_mzbach)
}

summary(discolon[[1]])

dfs<-discolon[[1]]
scores<-c("BMI_p5", "PI_mdson_sc",  "PI_mc_sc",  
          "PI_hmn_sc", "PI_do_sc",  "PI_mzbach_sc")

outcomes <- c('OS', 'DOC', 'CSS')


blood_cox<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxfull(data = dfs, entrytime = 'time_blood', score = score, outcome = outcome)
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

library(readxl)
multicox_subgroup <- read_excel("methy_BMI/results/multicox_subgroup.xlsx")


result_wide$Subgroup<-'Distal colon'
result_wide<-result_wide[, c(1, 8, 6, 4, 2)]
colnames(result_wide)<-colnames(multicox_subgroup)

multicox_subgroup<-rbind(multicox_subgroup, result_wide)
  

write_xlsx(multicox_subgroup, 
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/multicox_subgroup.xlsx')

agel50<-subset(blood_clin, Age_diag<50)  ## 98 
summary(agel50)


summary(blood_clin)








