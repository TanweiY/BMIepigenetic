source("C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)
library(survival)
library(writexl)

load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
## only use one dataset, because the percentages of missing values are not too high
df<-blood_all[[9]]
summary(df)
dput(names(df))
df$time_blood<-df$time_blood/30
summary(df$time_blood)

############## collect <3 months ##############
df_e<-subset(df, time_blood<3)

df_e<-within.data.frame(df_e, {
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
  
  BMI_cat<- relevel(BMI_cat, ref = '20-25')
  BMI5_14er_cat<- relevel(BMI5_14er_cat, ref = '20-25')
  meanBMI_tdiagcat<- relevel(meanBMI_tdiagcat, ref = '20-25')
  meanBMI_t514cat<- relevel(meanBMI_t514cat, ref = '20-25')
  
})

summary(df_e)

# Define the scores and BMI measures
cont_var<-c("PI_mc", "PI_hmn", "PI_do", "PI_mdson", "PI_mzbach")

df_e[cont_var]<-scale(df_e[cont_var])


scores <- c("PI_mdson", "PI_mdsonqt", "PI_mc", "PI_mcqt", "PI_hmn","PI_hmnqt", "PI_do","PI_doqt",  "PI_mzbach", "PI_mzbachqt")

outcomes<-c('OS', 'DOC', 'CSS', 'TTR')

blood_cox<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxfull(data = df_e, entrytime = 'time_blood', score = score, outcome = outcome)
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

result_wide$Group<-'Blood collected < 3 months'
result_wide$Sample<-'1233'
write_xlsx(result_wide,
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/cox_btimel3m.xlsx')

############## collect >3 months ##############

df_e<-subset(df, time_blood>=3) 

df_e<-within.data.frame(df_e, {
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
  
  BMI_cat<- relevel(BMI_cat, ref = '20-25')
  BMI5_14er_cat<- relevel(BMI5_14er_cat, ref = '20-25')
  meanBMI_tdiagcat<- relevel(meanBMI_tdiagcat, ref = '20-25')
  meanBMI_t514cat<- relevel(meanBMI_t514cat, ref = '20-25')
  
})

summary(df_e)

# Define the scores and BMI measures
cont_var<-c("PI_mc", "PI_hmn", "PI_do", "PI_mdson", "PI_mzbach")

df_e[cont_var]<-scale(df_e[cont_var])


scores <- c("PI_mdson", "PI_mdsonqt", "PI_mc", "PI_mcqt", "PI_hmn","PI_hmnqt", "PI_do","PI_doqt",  "PI_mzbach", "PI_mzbachqt")

outcomes<-c('OS', 'DOC', 'CSS', 'TTR')

blood_cox<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxfull(data = df_e, entrytime = 'time_blood', score = score, outcome = outcome)
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

result_wide$Group<-'Blood collected > 3 months'
result_wide$Sample<-'893'
write_xlsx(result_wide,
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/cox_btimem3m.xlsx')


