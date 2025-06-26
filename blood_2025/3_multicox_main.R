library(cgwtools)
library(dplyr)
library(writexl)

############## prepare the dataset first #############
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
summary(blood_all[[1]])
for (i in 1:20) {
  
  #### continous variables 
  blood_all[[i]]$BMI_p5<-blood_all[[i]]$BMI/5
  blood_all[[i]]$BMIt514_p5<-blood_all[[i]]$BMI_5_14earlier/5
  blood_all[[i]]$meanBMItd_p5<-blood_all[[i]]$meanBMI_tdiag/5
  blood_all[[i]]$meanBMItd_514<-blood_all[[i]]$meanBMI_t514/5
  
  blood_all[[i]]$PI_mdson_sc<-scale(blood_all[[i]]$PI_mdson)
  blood_all[[i]]$PI_mc_sc<-scale(blood_all[[i]]$PI_mc)
  blood_all[[i]]$PI_hmn_sc<-scale(blood_all[[i]]$PI_hmn)
  blood_all[[i]]$PI_do_sc<-scale(blood_all[[i]]$PI_do)
  blood_all[[i]]$PI_mzbach_sc<-scale(blood_all[[i]]$PI_mzbach)

  ### categorize variables
  blood_all[[i]]$BMI_cat<-cut(blood_all[[i]]$BMI, 
                                  breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  blood_all[[i]]$BMI_cat<-relevel(blood_all[[i]]$BMI_cat,
                                      ref = '20-25')
  
  blood_all[[i]]$BMI5_14er_cat<- cut(blood_all[[i]]$BMI_5_14earlier, 
                                         breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  
  blood_all[[i]]$BMI5_14er_cat<-relevel(blood_all[[i]]$BMI5_14er_cat,
                                            ref = '20-25')
  
  blood_all[[i]]$meanBMI_tdiagcat<- cut(blood_all[[i]]$meanBMI_tdiag, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  blood_all[[i]]$meanBMI_tdiagcat<-relevel(blood_all[[i]]$meanBMI_tdiagcat,
                                               ref = '20-25')
  
  blood_all[[i]]$meanBMI_t514cat<- cut(blood_all[[i]]$meanBMI_t514, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  blood_all[[i]]$meanBMI_t514cat<-relevel(blood_all[[i]]$meanBMI_t514cat,
                                              ref = '20-25')
  
  blood_all[[i]]$PI_mdsonqt<-factor(ntile(blood_all[[i]]$PI_mdson, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  blood_all[[i]]$PI_mcqt<-factor(ntile(blood_all[[i]]$PI_mc, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  blood_all[[i]]$PI_hmnqt<-factor(ntile(blood_all[[i]]$PI_hmn, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  blood_all[[i]]$PI_doqt<-factor(ntile(blood_all[[i]]$PI_do, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
  blood_all[[i]]$PI_mzbachqt<-factor(ntile(blood_all[[i]]$PI_mzbach, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"))
}

summary(blood_all[[1]])
save(blood_all, 
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")

#### start the multicox regression for various measurements ############
source('C:/Users/T532N/Desktop/BMI_m/methy_BMI/code/function/Cox_imputed_function.R')
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
dput(names(blood_all[[1]]))
dfs<-blood_all
scores<-c("BMI_p5","BMI_cat", "BMIt514_p5", "BMI5_14er_cat",
          "meanBMItd_p5", "meanBMI_tdiagcat", "meanBMItd_514","meanBMI_t514cat",
          "PI_mdson_sc", "PI_mdsonqt",  "PI_mc_sc", "PI_mcqt",  "PI_hmn_sc", "PI_hmnqt",
          "PI_do_sc", "PI_doqt", "PI_mzbach_sc",   "PI_mzbachqt")

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
           'C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/Multicox_allsample.xlsx')





























