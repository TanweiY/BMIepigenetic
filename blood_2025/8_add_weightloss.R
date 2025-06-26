library(cgwtools)
library(dplyr)
library(writexl)

############## prepare the dataset first #############
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")

summary(blood_all[[1]])
for (i in 1:20) {
  blood_all[[i]]$height<-blood_all[[i]]$height/100
  blood_all[[i]]$weight_diag<- blood_all[[i]]$BMI * blood_all[[i]]$height^2
  blood_all[[i]]$weight_514yr<-blood_all[[i]]$BMI_5_14earlier * blood_all[[i]]$height^2
  blood_all[[i]]$weight_514change<- blood_all[[i]]$weight_diag - blood_all[[i]]$weight_514yr
  }

summary(blood_all[[1]])

# 1183
for (i in 1:20) {
  blood_all[[i]]$weight_514changecat<-cut(blood_all[[i]]$weight_514change, 
                                          breaks = c(-Inf, -5, -2, Inf), labels = c('> 5kg','2-5 kg','no weight loss'))
  
  blood_all[[i]]$weight_514changecat<-relevel(blood_all[[i]]$weight_514changecat, ref = 'no weight loss')
}

summary(blood_all[[i]]$weight_514changecat)

save(blood_all, 
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")

#### start the multicox regression for weigt loss ############
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/Cox_imputed_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
dput(names(blood_all[[1]]))
dfs<-blood_all
scores<-c("weight_514changecat")
outcomes <- c('OS', 'DOC', 'CSS')

multicox_results<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_fullmulticox_impu(dfs,  entrytime = 'time_blood', score = score, outcome = outcome)
    multicox_results<- rbind(multicox_results, result)
  }
}

multicox_results<-multicox_results[order(multicox_results$Outcome), ]
multicox_results<-multicox_results[, c(1, 6, 7)]
multicox_results<-reshape(multicox_results, 
                          timevar = "Outcome",
                          idvar = "Level",
                          direction = "wide",
                          v.names = c("HR_95CI"))

write_xlsx(multicox_results, 
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/coxweightloss.xlsx')


library(cgwtools)
library(dplyr)
library(writexl)

############## prepare the dataset first #############
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")

summary(blood_all[[1]])
for (i in 1:20) {
  blood_all[[i]]$height<-blood_all[[i]]$height/100
  blood_all[[i]]$weight_diag<- blood_all[[i]]$BMI * blood_all[[i]]$height^2
  blood_all[[i]]$weight_514yr<-blood_all[[i]]$BMI_5_14earlier * blood_all[[i]]$height^2
  blood_all[[i]]$weight_514change<- blood_all[[i]]$weight_diag - blood_all[[i]]$weight_514yr
}

summary(blood_all[[1]])
for (i in 1:20) {
  blood_all[[i]]<-subset(blood_all[[i]], weight_514change < -2)
}

# 943
save(blood_all, 
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/weightlossev20.RData")

#### start the multicox regression for various measurements ############
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/Cox_imputed_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/noweightchangev20.RData")
dput(names(blood_all[[1]]))
dfs<-blood_all
scores<-c("BMI_p5","BMI_cat", 
          "PI_mdson_sc", "PI_mdsonqt",  "PI_mc_sc", "PI_mcqt",  "PI_hmn_sc", "PI_hmnqt",
          "PI_do_sc", "PI_doqt", "PI_mzbach_sc",   "PI_mzbachqt")

outcomes <- c('OS', 'DOC', 'CSS')

weightlosss<-data.frame()
for (score in scores) {
  for (outcome in outcomes) {
    result <- run_fullmulticox_impu(dfs,  entrytime = 'time_blood', score = score, outcome = outcome)
    weightlosss<- rbind(weightlosss, result)
  }
}

weightlosss<-weightlosss[order(weightlosss$Outcome), ]
weightlosss<-weightlosss[, c(1, 6, 7)]
weightlosss<-reshape(weightlosss, 
                     timevar = "Outcome",
                     idvar = "Level",
                     direction = "wide",
                     v.names = c("HR_95CI"))

write_xlsx(weightlosss, 
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/weightloss.xlsx')


noweightloss$subgroup<-'no weight loss (weight change >-2 kg)'
weightlosss$subgroup<-'weight loss (weight change <-2 kg)'

weightchangesubgroup<-



















































