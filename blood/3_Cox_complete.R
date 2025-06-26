source("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_clin.RData")

### stage I-III
blood_clin<-subset(blood_clin, TNM_adj !='IV') ## 1909

summary(blood_clin)

dput(names(blood_clin))
var<-c("Age_diag", "Sex" , "TNM_adj" ,"BMI", "high_blood_pressure", "NSAIDS",
       "smoking", "cardiovad", "statins", "chemradther", "crc2sites", "physical_activity_lifeav",
       "evercol", "hormone_replace")

var<-c("Age_diag", "Sex" , "TNM_adj")

## complete case
blood_clin_complet<-blood_clin[complete.cases(blood_clin[var]), ]
## 1830 --> only 4% missing

blood_clin<-within.data.frame(blood_clin, {
  PI_epitob_md<-factor(ifelse(PI_epitob < median(PI_epitob), 'Low', 'High'), levels = c('Low', 'High'))
  PI_mc_md<-factor(ifelse(PI_mc < median(PI_mc), 'Low', 'High'), levels = c('Low', 'High'))
  PI_Liu_md<-factor(ifelse(PI_liu < median(PI_liu), 'Low', 'High'), levels = c('Low', 'High'))
  PI_epitob_tert<-factor(ifelse(PI_epitob < quantile(PI_epitob, 0.33), 'T1', 
                                ifelse(PI_epitob > quantile(PI_epitob, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  PI_mc_tert<-factor(ifelse(PI_mc < quantile(PI_mc, 0.33), 'T1', 
                                ifelse(PI_mc > quantile(PI_mc, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  PI_Liu_tert<-factor(ifelse(PI_liu < quantile(PI_liu, 0.33), 'T1', 
                            ifelse(PI_liu > quantile(PI_liu, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  liftimdrink_cat<-relevel(liftimdrink_cat, ref = 'Light drinkers')
  
  recentdrink_cat<-relevel(recentdrink_cat, ref = 'Light drinkers')

  
  })

dput(names(blood_clin))

score <- c("PI_epitob", "PI_mc", "PI_liu")

blood_clin[score]<-scale(blood_clin[score])

scores <- c("liftimdrink_cat", "recentdrink_cat", "PI_epitob", "PI_mc", "PI_liu",
            "PI_Liu_md", "PI_mc_md", "PI_epitob_md", 
            "PI_epitob_tert", "PI_Liu_tert", "PI_mc_tert")

outcomes<-c('OS', 'DOC', 'CSS')

blood_cox<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxfull(data = blood_clin, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox<- rbind(blood_cox, result)
  }
}
blood_cox$dataset<-'complete cases'
blood_cox$adjustment<-'fully adjusted'
blood_cox$stage<-'I-III'

blood_cox_part<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxpart(data = blood_clin, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox_part<- rbind(blood_cox_part, result)
  }
}
blood_cox_part$dataset<-'complete cases'
blood_cox_part$adjustment<-'age, sex, stage'
blood_cox_part$stage<-'I-III'

blood_completecox<-bind_rows(blood_cox, blood_cox_part)

blood_completecox$`aHR(95% CI)`<-gsub("\\[", "(", blood_completecox$`aHR(95% CI)`)
blood_completecox$`aHR(95% CI)`<-gsub("\\]", ")", blood_completecox$`aHR(95% CI)`)

write_xlsx(blood_completecox, 
       '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/blood_completecox.xlsx')



