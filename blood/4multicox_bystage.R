source("C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)
library(survival)

load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
## only use one dataset, because the percentages of missing values are not too high
df<-blood_all[[9]]
summary(df)
dput(names(df))

############## stage I-III ##############
df<-blood_all[[9]]
summary(df)
dput(names(df))
summary(df$TNM_adj)
df<-subset(df, TNM_adj =='I'|TNM_adj =='II'|TNM_adj =='III')
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
  
  BMI_cat<- relevel(BMI_cat, ref = '20-25')
  BMI5_14er_cat<- relevel(BMI5_14er_cat, ref = '20-25')
  meanBMI_tdiagcat<- relevel(meanBMI_tdiagcat, ref = '20-25')
  meanBMI_t514cat<- relevel(meanBMI_t514cat, ref = '20-25')
  
})

summary(df)

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

result_wide$Group<-'StageI_III'
result_wide$Sample<-'1832'
s1_3<-result_wide

write_xlsx(s1_3,
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/cox_stage1_3.xlsx')

###### only stage IV ####
source("C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)
library(survival)

load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
## only use one dataset, because the percentages of missing values are not too high
df<-blood_all[[9]]
summary(df)
dput(names(df))

############## stage IV ##############
df<-blood_all[[9]]
summary(df)
dput(names(df))
summary(df$TNM_adj)
df<-subset(df, TNM_adj =='IV')
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
  
  BMI_cat<- relevel(BMI_cat, ref = '20-25')
  BMI5_14er_cat<- relevel(BMI5_14er_cat, ref = '20-25')
  meanBMI_tdiagcat<- relevel(meanBMI_tdiagcat, ref = '20-25')
  meanBMI_t514cat<- relevel(meanBMI_t514cat, ref = '20-25')
  
})

summary(df)

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

result_wide$Group<-'StageIV'
result_wide$Sample<-'294'
s4<-result_wide

write_xlsx(s4,
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/cox_stage4.xlsx')

