source("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/code/Cox_imputed_function.R")
library(cgwtools)
library(dplyr)
library(writexl)

####################### run the imputed dataset ####################
## full case
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_imputatedpro.RData")
dput(names(blood_imputatedpro[[1]]))
## run the same function for stage I-III and define the median, tertile cut-off
for (i in 1:20){
  blood_imputatedpro[[i]]<-subset(blood_imputatedpro[[i]], 
                                  TNM_adj != 'IV') ## 1918

  
  blood_imputatedpro[[i]]$PI_epitob_md<-factor(ifelse(blood_imputatedpro[[i]]$PI_epitob < median(blood_imputatedpro[[i]]$PI_epitob), 
                                                      'Low', 'High'), levels = c('Low', 'High'))

  blood_imputatedpro[[i]]$PI_Liu_md<-factor(ifelse(blood_imputatedpro[[i]]$PI_liu < median(blood_imputatedpro[[1]]$PI_liu), 
                                                      'Low', 'High'), levels = c('Low', 'High'))
  
  blood_imputatedpro[[i]]$PI_mc_md<-factor(ifelse(blood_imputatedpro[[i]]$PI_mc < median(blood_imputatedpro[[i]]$PI_mc), 
                                                   'Low', 'High'), levels = c('Low', 'High'))
  
  blood_imputatedpro[[i]]$PI_epitob_tert<-factor(ifelse(blood_imputatedpro[[i]]$PI_epitob < quantile(blood_imputatedpro[[i]]$PI_epitob, 0.33), 'T1', 
                                ifelse(blood_imputatedpro[[i]]$PI_epitob > quantile(blood_imputatedpro[[i]]$PI_epitob, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  blood_imputatedpro[[i]]$PI_Liu_tert<-factor(ifelse(blood_imputatedpro[[i]]$PI_liu < quantile(blood_imputatedpro[[i]]$PI_liu, 0.33), 'T1', 
                                                        ifelse(blood_imputatedpro[[i]]$PI_liu > quantile(blood_imputatedpro[[i]]$PI_liu, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  blood_imputatedpro[[i]]$PI_mc_tert<-factor(ifelse(blood_imputatedpro[[i]]$PI_mc < quantile(blood_imputatedpro[[i]]$PI_mc, 0.33), 'T1', 
                                                     ifelse(blood_imputatedpro[[i]]$PI_mc > quantile(blood_imputatedpro[[i]]$PI_mc, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  blood_imputatedpro[[i]]$liftimdrink_cat<-relevel(blood_imputatedpro[[i]]$liftimdrink_cat, ref = 'Light drinkers')
  
  blood_imputatedpro[[i]]$recentdrink_cat<-relevel(blood_imputatedpro[[i]]$recentdrink_cat, ref = 'Light drinkers')
  
  
}

dput(names(blood_imputatedpro[[1]]))
dfs <- blood_imputatedpro # List of 20 dataframes
scores <- c("liftimdrink_cat", "recentdrink_cat", "PI_epitob", "PI_mc", "PI_liu",
             "PI_mc_md", "PI_epitob_md", "PI_Liu_md",
            "PI_mc_tert", "PI_epitob_tert", "PI_Liu_tert")

outcomes <- c('OS', 'DOC', 'CSS')

# Run the function
blood_cox_part<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_partmulticox_impu(dfs, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox_part<- rbind(blood_cox_part, result)
  }
}

blood_cox_part<-blood_cox_part[order(blood_cox_part$Level, blood_cox_part$Outcome), ]  
blood_cox_part<-blood_cox_part[, c(1,6,7)]
blood_cox_part$dataset<-'20 imputed datasets'
blood_cox_part$adjustment<-'age, sex, stage'
blood_cox_part$stage<-'I-III'

blood_cox_full<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_fullmulticox_impu(dfs, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox_full<- rbind(blood_cox_full, result)
  }
}

blood_cox_full<-blood_cox_full[order(blood_cox_full$Level, blood_cox_full$Outcome), ]  
blood_cox_full<-blood_cox_full[, c(1,6,7)]
blood_cox_full$dataset<-'20 imputed datasets'
blood_cox_full$adjustment<-'fully adjusted'
blood_cox_full$stage<-'I-III'

bloodimpucox13<-bind_rows(blood_cox_part,
                          blood_cox_full)

write_xlsx(bloodimpucox13, 
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/bloodimpucox13.xlsx')

## run the same function for stage IV
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_imputatedpro.RData")

## run the same function for stage IV and define the median, tertile cut-off
for (i in 1:20){
  blood_imputatedpro[[i]]<-subset(blood_imputatedpro[[i]], 
                                  TNM_adj == 'IV') ## 1918
  
  blood_imputatedpro[[i]]$PI_epitob_md<-factor(ifelse(blood_imputatedpro[[i]]$PI_epitob < median(blood_imputatedpro[[i]]$PI_epitob), 
                                                      'Low', 'High'), levels = c('Low', 'High'))
  
  blood_imputatedpro[[i]]$PI_Liu_md<-factor(ifelse(blood_imputatedpro[[i]]$PI_liu < median(blood_imputatedpro[[1]]$PI_liu), 
                                                   'Low', 'High'), levels = c('Low', 'High'))
  
  blood_imputatedpro[[i]]$PI_mc_md<-factor(ifelse(blood_imputatedpro[[i]]$PI_mc < median(blood_imputatedpro[[i]]$PI_mc), 
                                                  'Low', 'High'), levels = c('Low', 'High'))
  
  blood_imputatedpro[[i]]$PI_epitob_tert<-factor(ifelse(blood_imputatedpro[[i]]$PI_epitob < quantile(blood_imputatedpro[[i]]$PI_epitob, 0.33), 'T1', 
                                                        ifelse(blood_imputatedpro[[i]]$PI_epitob > quantile(blood_imputatedpro[[i]]$PI_epitob, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  blood_imputatedpro[[i]]$PI_Liu_tert<-factor(ifelse(blood_imputatedpro[[i]]$PI_liu < quantile(blood_imputatedpro[[i]]$PI_liu, 0.33), 'T1', 
                                                     ifelse(blood_imputatedpro[[i]]$PI_liu > quantile(blood_imputatedpro[[i]]$PI_liu, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  blood_imputatedpro[[i]]$PI_mc_tert<-factor(ifelse(blood_imputatedpro[[i]]$PI_mc < quantile(blood_imputatedpro[[i]]$PI_mc, 0.33), 'T1', 
                                                    ifelse(blood_imputatedpro[[i]]$PI_mc > quantile(blood_imputatedpro[[i]]$PI_mc, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
}

dput(names(blood_imputatedpro[[1]]))
dfs <- blood_imputatedpro # List of 20 dataframes
scores <- c("liftimdrink_cat", "recentdrink_cat", "PI_epitob", "PI_mc", "PI_liu",
            "PI_epitob_md", "PI_Liu_md", "PI_mc_md",
            "PI_epitob_tert", "PI_Liu_tert", "PI_mc_tert")

outcomes <- c('OS', 'DOC', 'CSS')

# Run the function
blood_cox_part<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_partmulticox_impu(dfs, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox_part<- rbind(blood_cox_part, result)
  }
}

blood_cox_part<-blood_cox_part[order(blood_cox_part$Level, blood_cox_part$Outcome), ]  
blood_cox_part<-blood_cox_part[, c(1,6,7)]
blood_cox_part$dataset<-'20 imputed datasets'
blood_cox_part$adjustment<-'age, sex, stage'
blood_cox_part$stage<-'IV'

blood_cox_full<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_fullmulticox_impu(dfs, entrytime = 'time_blood', score = score, outcome = outcome)
    blood_cox_full<- rbind(blood_cox_full, result)
  }
}

blood_cox_full<-blood_cox_full[order(blood_cox_full$Level, blood_cox_full$Outcome), ]  
blood_cox_full<-blood_cox_full[, c(1,6,7)]
blood_cox_full$dataset<-'20 imputed datasets'
blood_cox_full$adjustment<-'fully adjusted'
blood_cox_full$stage<-'IV'

bloodimpucox4<-bind_rows(blood_cox_part,
                         blood_cox_full)

write_xlsx(bloodimpucox4, 
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/bloodimpucox4.xlsx')





  