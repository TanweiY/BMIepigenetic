source("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/code/Cox_complete_function.R")
library(cgwtools)
library(dplyr)
library(tableone)
library(survival)

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")
summary(tumor_clin$liftimdrink_cat)
summary(tumor_clin$recentdrink_cat)

### stage I-III
tumor_clin<-subset(tumor_clin, TNM_adj !='IV') ## 1840

summary(tumor_clin)

dput(names(tumor_clin))
var<-c("Age_diag", "Sex" , "TNM_adj" ,"BMI", "high_blood_pressure", "NSAIDS",
       "smoking", "cardiovad", "statins", "chemradther", "crc2sites", "physical_activity_lifeav",
       "evercol", "hormone_replace")

var<-c("Age_diag", "Sex" , "TNM_adj")
  

## complete case
tumor_clin_complet<-tumor_clin[complete.cases(tumor_clin[var]), ] # 1840


## 1830 --> only 4% missing

tumor_clin<-within.data.frame(tumor_clin, {
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

summary(tumor_clin$liftimdrink_cat)
summary(tumor_clin$recentdrink_cat)


dput(names(tumor_clin))

score <- c("PI_epitob", "PI_mc", "PI_liu")

tumor_clin[score]<-scale(tumor_clin[score])

scores <- c("liftimdrink_cat", "recentdrink_cat", "PI_epitob", "PI_mc", "PI_liu",
            "PI_Liu_md", "PI_mc_md", "PI_epitob_md", 
            "PI_epitob_tert", "PI_Liu_tert", "PI_mc_tert")

outcomes<-c('OS', 'DOC', 'CSS')

tumor_cox<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxfull(data = tumor_clin, entrytime = 'late_days', score = score, outcome = outcome)
    tumor_cox<- rbind(tumor_cox, result)
  }
}
tumor_cox$dataset<-'complete cases'
tumor_cox$adjustment<-'fully adjusted'
tumor_cox$stage<-'I-III'

tumor_cox_part<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- multicoxpart(data = tumor_clin, entrytime = 'late_days', score = score, outcome = outcome)
    tumor_cox_part<- rbind(tumor_cox_part, result)
  }
}
tumor_cox_part$dataset<-'complete cases'
tumor_cox_part$adjustment<-'age, sex, stage'
tumor_cox_part$stage<-'I-III'

tumor_completecox<-bind_rows(tumor_cox, tumor_cox_part)

tumor_completecox$`aHR(95% CI)`<-gsub("\\[", "(", tumor_completecox$`aHR(95% CI)`)
tumor_completecox$`aHR(95% CI)`<-gsub("\\]", ")", tumor_completecox$`aHR(95% CI)`)

write_xlsx(tumor_completecox, 
       '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/tumor_completecox.xlsx')



