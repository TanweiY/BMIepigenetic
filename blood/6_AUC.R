library(dplyr)
library(writexl)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_clin.RData")
summary(blood_clin)
blood_clin<-within.data.frame(blood_clin, {
  recentdrink_binary<-factor(ifelse(lastyralcohol_day ==0, 'non-drinkers', 
                                    ifelse(is.na(lastyralcohol_day), NA, 'drinkers')), 
                             levels = c('non-drinkers', 'drinkers'))
  
  lifetimedrink_binary<-factor(ifelse(lifetimealcohol_day ==0, 'non-drinkers', 
                                    ifelse(is.na(lifetimealcohol_day), NA, 'drinkers')), 
                             levels = c('non-drinkers', 'drinkers'))
  
})

save(blood_clin,
     file ="/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_clin.RData")

## assess the AUC
library(pROC)
get_aucci<-function(dataset, outcome_col, score_col) {
  outcome<-dataset[[outcome_col]]
  score<-dataset[[score_col]]
  r<-roc(outcome, score)
  ci<-ci.auc(r, method = c('bootstrap'), conf.level = 0.95)
  result<-data.frame(Score = score_col, 
                     Outcome = outcome_col,
                     AUCCI = paste0(round(r$auc, 2), ' (', round(ci[1], 2), ', ',  round(ci[2], 2), ')'))
  return (result)
  
}

source('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/code/AUC_CI_function.R')
outcomes<-c('recentdrink_binary', 'lifetimedrink_binary')
scores <- c("PI_epitob", "PI_mc", "PI_liu")

AUC_result<-data.frame()
for (score in scores) {
  for (outcome in outcomes){
    result<-get_aucci(dataset = blood_clin, outcome_col = outcome, score_col = score)
    AUC_result<-rbind(AUC_result, result)
    
  }
}

AUC_result$source<-'blood'
## run the function on all scores 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")
summary(tumor_clin)
tumor_clin<-within.data.frame(tumor_clin, {
  recentdrink_binary<-factor(ifelse(lastyralcohol_day ==0, 'non-drinkers', 
                                    ifelse(is.na(lastyralcohol_day), NA, 'drinkers')), 
                             levels = c('non-drinkers', 'drinkers'))
  
  lifetimedrink_binary<-factor(ifelse(lifetimealcohol_day ==0, 'non-drinkers', 
                                      ifelse(is.na(lifetimealcohol_day), NA, 'drinkers')), 
                               levels = c('non-drinkers', 'drinkers'))
  
})

save(tumor_clin,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata"
     
     
     outcomes<-c('recentdrink_binary', 'lifetimedrink_binary')
     scores <- c("PI_epitob", "PI_mc", "PI_liu")
     
     AUC_tumor<-data.frame()
     for (score in scores) {
       for (outcome in outcomes){
         result<-get_aucci(dataset = tumor_clin, outcome_col = outcome, score_col = score)
         AUC_tumor<-rbind(AUC_tumor, result)
         
       }
     }
AUC_tumor$source<-'tumor'
AUC_all<-bind_rows(AUC_result, AUC_tumor)
write_xlsx(AUC_all,
      '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/AUC_all.xlsx')

     
