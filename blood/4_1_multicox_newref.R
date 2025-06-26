source("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)
library(survival)

load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_progdf.RData")
summary(df)

df<-within.data.frame(df, {
  BMI_cat<- relevel(BMI_cat, ref = '20-25')
  BMI5_14er_cat<- relevel(BMI5_14er_cat, ref = '20-25')
  meanBMI_tdiagcat<- relevel(meanBMI_tdiagcat, ref = '20-25')
  meanBMI_t514cat<- relevel(meanBMI_t514cat, ref = '20-25')
  
  WYO_tdqt<-relevel(WYO_tdqt, ref = 'Q2')
  WYO_tr5yrqt<-relevel(WYO_tr5yrqt, ref = 'Q2')
  PI_mdsonqt<-relevel(PI_mdsonqt, ref = 'Q2')
  PI_mcqt<-relevel(PI_mcqt, ref = 'Q2')
  PI_hmnqt<-relevel(PI_hmnqt, ref = 'Q2')
  PI_doqt<-relevel(PI_doqt, ref = 'Q2')
  PI_mzbachqt<-relevel(PI_mzbachqt, ref = 'Q2')
})

summary(df)

save(df,
     file = '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_progrelevel.RData')

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

write_xlsx(result_wide,
           '/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/results/cox_results_relevel.xlsx')







