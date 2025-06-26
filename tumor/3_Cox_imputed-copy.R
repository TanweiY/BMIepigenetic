source("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/code/Cox_imputed_function.R")
library(cgwtools)
library(dplyr)
library(writexl)

####################### run the imputed dataset ####################
## full case
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_imputatedpro.RData")
dput(names(tumor_imputatedpro[[1]]$la))
## run the same function for stage I-III and define the median, tertile cut-off
for (i in 1:20){
  tumor_imputatedpro[[i]]<-subset(tumor_imputatedpro[[i]], 
                                  TNM_adj != 'IV') ## 1946
  
  
  tumor_imputatedpro[[i]]$PI_epitob_md<-factor(ifelse(tumor_imputatedpro[[i]]$PI_epitob < median(tumor_imputatedpro[[i]]$PI_epitob), 
                                                      'Low', 'High'), levels = c('Low', 'High'))
  
  tumor_imputatedpro[[i]]$PI_Liu_md<-factor(ifelse(tumor_imputatedpro[[i]]$PI_liu < median(tumor_imputatedpro[[1]]$PI_liu), 
                                                   'Low', 'High'), levels = c('Low', 'High'))
  
  tumor_imputatedpro[[i]]$PI_mc_md<-factor(ifelse(tumor_imputatedpro[[i]]$PI_mc < median(tumor_imputatedpro[[i]]$PI_mc), 
                                                  'Low', 'High'), levels = c('Low', 'High'))
  
  tumor_imputatedpro[[i]]$PI_epitob_tert<-factor(ifelse(tumor_imputatedpro[[i]]$PI_epitob < quantile(tumor_imputatedpro[[i]]$PI_epitob, 0.33), 'T1', 
                                                        ifelse(tumor_imputatedpro[[i]]$PI_epitob > quantile(tumor_imputatedpro[[i]]$PI_epitob, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  tumor_imputatedpro[[i]]$PI_Liu_tert<-factor(ifelse(tumor_imputatedpro[[i]]$PI_liu < quantile(tumor_imputatedpro[[i]]$PI_liu, 0.33), 'T1', 
                                                     ifelse(tumor_imputatedpro[[i]]$PI_liu > quantile(tumor_imputatedpro[[i]]$PI_liu, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  tumor_imputatedpro[[i]]$PI_mc_tert<-factor(ifelse(tumor_imputatedpro[[i]]$PI_mc < quantile(tumor_imputatedpro[[i]]$PI_mc, 0.33), 'T1', 
                                                    ifelse(tumor_imputatedpro[[i]]$PI_mc > quantile(tumor_imputatedpro[[i]]$PI_mc, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  
  tumor_imputatedpro[[i]]$liftimdrink_cat<-relevel(tumor_imputatedpro[[i]]$liftimdrink_cat, ref = 'Light drinkers')
  
  tumor_imputatedpro[[i]]$recentdrink_cat<-relevel(tumor_imputatedpro[[i]]$recentdrink_cat, ref = 'Light drinkers')
}

summary(tumor_imputatedpro[[1]]$liftimdrink_cat)


dput(names(tumor_imputatedpro[[1]]))
dfs <- tumor_imputatedpro # List of 20 dataframes
scores <- c("liftimdrink_cat", "recentdrink_cat", "PI_epitob", "PI_mc", "PI_liu",
            "PI_mc_md", "PI_epitob_md", "PI_Liu_md",
            "PI_mc_tert", "PI_epitob_tert", "PI_Liu_tert")

outcomes <- c('OS', 'DOC', 'CSS')

# Run the function
tumor_cox_part<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_partmulticox_impu(dfs, entrytime = 'late_days', score = score, outcome = outcome)
    tumor_cox_part<- rbind(tumor_cox_part, result)
  }
}

tumor_cox_part<-tumor_cox_part[order(tumor_cox_part$Level, tumor_cox_part$Outcome), ]  
tumor_cox_part<-tumor_cox_part[, c(1,6,7)]
tumor_cox_part$dataset<-'20 imputed datasets'
tumor_cox_part$adjustment<-'age, sex, stage'
tumor_cox_part$stage<-'I-III'

tumor_cox_full<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_fullmulticox_impu(dfs, entrytime = 'late_days', score = score, outcome = outcome)
    tumor_cox_full<- rbind(tumor_cox_full, result)
  }
}

tumor_cox_full<-tumor_cox_full[order(tumor_cox_full$Level, tumor_cox_full$Outcome), ]  
tumor_cox_full<-tumor_cox_full[, c(1,6,7)]
tumor_cox_full$dataset<-'20 imputed datasets'
tumor_cox_full$adjustment<-'fully adjusted'
tumor_cox_full$stage<-'I-III'

tumorimpucox13<-bind_rows(tumor_cox_part,
                          tumor_cox_full)

write_xlsx(tumorimpucox13, 
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/tumorimpucox13.xlsx')

## run the same function for stage IV
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_imputatedpro.RData")

## run the same function for stage IV and define the median, tertile cut-off
for (i in 1:20){
  tumor_imputatedpro[[i]]<-subset(tumor_imputatedpro[[i]], 
                                  TNM_adj == 'IV') ## 1918
  
  tumor_imputatedpro[[i]]$PI_epitob_md<-factor(ifelse(tumor_imputatedpro[[i]]$PI_epitob < median(tumor_imputatedpro[[i]]$PI_epitob), 
                                                      'Low', 'High'), levels = c('Low', 'High'))
  
  tumor_imputatedpro[[i]]$PI_Liu_md<-factor(ifelse(tumor_imputatedpro[[i]]$PI_liu < median(tumor_imputatedpro[[1]]$PI_liu), 
                                                   'Low', 'High'), levels = c('Low', 'High'))
  
  tumor_imputatedpro[[i]]$PI_mc_md<-factor(ifelse(tumor_imputatedpro[[i]]$PI_mc < median(tumor_imputatedpro[[i]]$PI_mc), 
                                                  'Low', 'High'), levels = c('Low', 'High'))
  
  tumor_imputatedpro[[i]]$PI_epitob_tert<-factor(ifelse(tumor_imputatedpro[[i]]$PI_epitob < quantile(tumor_imputatedpro[[i]]$PI_epitob, 0.33), 'T1', 
                                                        ifelse(tumor_imputatedpro[[i]]$PI_epitob > quantile(tumor_imputatedpro[[i]]$PI_epitob, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  tumor_imputatedpro[[i]]$PI_Liu_tert<-factor(ifelse(tumor_imputatedpro[[i]]$PI_liu < quantile(tumor_imputatedpro[[i]]$PI_liu, 0.33), 'T1', 
                                                     ifelse(tumor_imputatedpro[[i]]$PI_liu > quantile(tumor_imputatedpro[[i]]$PI_liu, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
  tumor_imputatedpro[[i]]$PI_mc_tert<-factor(ifelse(tumor_imputatedpro[[i]]$PI_mc < quantile(tumor_imputatedpro[[i]]$PI_mc, 0.33), 'T1', 
                                                    ifelse(tumor_imputatedpro[[i]]$PI_mc > quantile(tumor_imputatedpro[[i]]$PI_mc, 0.67), 'T3', 'T2')), levels = c('T1', 'T2', 'T3'))
  
}

dput(names(tumor_imputatedpro[[1]]))
dfs <- tumor_imputatedpro # List of 20 dataframes
scores <- c("liftimdrink_cat", "recentdrink_cat", "PI_epitob", "PI_mc", "PI_liu",
            "PI_epitob_md", "PI_Liu_md", "PI_mc_md",
            "PI_epitob_tert", "PI_Liu_tert", "PI_mc_tert")

outcomes <- c('OS', 'DOC', 'CSS')

# Run the function
tumor_cox_part<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_partmulticox_impu(dfs, entrytime = 'late_days', score = score, outcome = outcome)
    tumor_cox_part<- rbind(tumor_cox_part, result)
  }
}

tumor_cox_part<-tumor_cox_part[order(tumor_cox_part$Level, tumor_cox_part$Outcome), ]  
tumor_cox_part<-tumor_cox_part[, c(1,6,7)]
tumor_cox_part$dataset<-'20 imputed datasets'
tumor_cox_part$adjustment<-'age, sex, stage'
tumor_cox_part$stage<-'IV'

tumor_cox_full<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_fullmulticox_impu(dfs, entrytime = 'late_days', score = score, outcome = outcome)
    tumor_cox_full<- rbind(tumor_cox_full, result)
  }
}

tumor_cox_full<-tumor_cox_full[order(tumor_cox_full$Level, tumor_cox_full$Outcome), ]  
tumor_cox_full<-tumor_cox_full[, c(1,6,7)]
tumor_cox_full$dataset<-'20 imputed datasets'
tumor_cox_full$adjustment<-'fully adjusted'
tumor_cox_full$stage<-'IV'

tumorimpucox4<-bind_rows(tumor_cox_part,
                         tumor_cox_full)

write_xlsx(tumorimpucox4, 
           '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/tumorimpucox4.xlsx')





