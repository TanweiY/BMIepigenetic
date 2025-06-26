library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(htmlTable)
library(rmarkdown)
library(devtools)
library(EpiSmokEr)
library(readxl)
library(cgwtools)
library(writexl)

library(minfi)
library(limma)
library(impute)
library(preprocessCore)
library(BiocParallel)
library(dplyr)

####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/bloodm_ILM.RData")
sum(is.na(dataset_ILM))
## perform imputation
beta_ILM<-impute.knn(dataset_ILM)$data
sum(is.na(beta_ILM))

scores <-  read_excel("methy_BMI/BMI_cpgs.xlsx")
cpg<-as.data.frame(rownames(beta_ILM)) # 5585

intersect<-Reduce(intersect, list(scores$CpGs, cpg$`rownames(beta_ILM)`)) # 4920
beta_ILM<-beta_ILM[intersect, ]

beta_ILM <- as.data.frame(beta_ILM)
beta_ILM$cpg<-rownames(beta_ILM)
beta_ILM[nrow(beta_ILM)+1, ]<-colnames(beta_ILM)
library(data.table)
beta_ILM<-transpose(beta_ILM)
colnames(beta_ILM)<-beta_ILM[nrow(beta_ILM), ]
names(beta_ILM)[ncol(beta_ILM)]<-'uid'
library(dplyr)
beta_ILM <- beta_ILM %>%
  select(uid, everything())
beta_ILM<-beta_ILM[-nrow(beta_ILM), ]
sum(is.na(beta_ILM))

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/blood_methylation/processed_data/methy_id.rdata")

beta_ILM<-merge(id, beta_ILM, by = 'uid')
beta_ILM[intersect]<-lapply(beta_ILM[intersect], function(x)as.numeric(x))

save(beta_ILM,
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/beta_ILM_blood.RData')


############## clinical variables preparation ###########
library(haven)
library(dplyr)
library(readxl)
baseline<-read_sas('C:/Users/T532N/Desktop/cumulative_BMI_prognosis/dachs_cohort/dachskomplett_categ_20191231.sas7bdat')

variables <- read_excel('C:/Users/T532N/Desktop/BMI_m/methy_BMI/variables.xlsx', 
                        sheet = 'baseline')

base<-subset(baseline, select = variables$var)
colnames(base)<-c(variables$Names)
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/beta_ILM_blood.RData")

##
base<-base[(base$lab_nummer %in% beta_ILM$tn_id), ] ## 2242
rm(beta_ILM)

## select follow up dataset
follow<-read_sas('C:/Users/T532N/Desktop/cumulative_BMI_prognosis/dachs_cohort/follow.sas7bdat')
variables <- read_excel("C:/Users/T532N/Desktop/BMI_m/methy_BMI/variables.xlsx", 
                        sheet = "follow")

f<-subset(follow, select = variables$var)
colnames(f)<-variables$Names

blood_clin<-merge(base, f, by = 'tn_id')

## exclude without lifetime alchohol consumption
summary(blood_clin$BMI) # 10 missing
summary(blood_clin$BMI_1_10earlier) # 41 missing
summary(blood_clin$BMI_5_14earlier) # 35 missing

blood_clin<-subset(blood_clin, !is.na(blood_clin$BMI))
summary(blood_clin)
blood_clin<-subset(blood_clin, !is.na(blood_clin$timeD)) ## 2231

## variable coding ##
summary(blood_clin)
colnames(blood_clin)
blood_clin<-within.data.frame(blood_clin, {
  time_blood<-as.numeric(difftime(Blutabnahme_Dat, indexdat,  units = "days"))
  BMI_cat<- cut(BMI, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  BMI1_10er_cat<- cut(BMI_1_10earlier, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  BMI5_14er_cat<- cut(BMI_5_14earlier, breaks = c(-Inf, 20, 25, 30, Inf), labels = c('<20','20-25','25-30', '≥30'))
  
  Schooling_years<- factor(Schooling_years, levels = c(1,2,3), labels = c('<9','9-10','>10'))
  
  Sex<-factor(Sex, levels = c(1,2), labels = c('Female', 'Male'))
  TNM_adj<-factor(TNM_adj, levels = c(1,2,3,4), labels = c('I', 'II', 'III', 'IV'), ordered = T)
  Location <- as.factor(ifelse(crc2sites=='rectum', 'Rectum', ifelse(crcprox == 1, 'Proximal colon', 'Distal colon')))
  smoking<-factor(smoking, levels = c(0, 1, 2), labels = c('Never', 'Former', 'Current'))
  statins<-factor(statins, levels = c(0,1), labels = c('No', 'Yes'))
  other_cancer<-factor(other_cancer, levels = c(0,1), labels = c('No', 'Yes'))
  high_blood_pressure<-factor(high_blood_pressure, levels = c(0,1), labels = c('No', 'Yes'))
  hyperlipidemia<-factor(hyperlipidemia, levels = c(0,1), labels = c('No', 'Yes'))
  NSAIDS<-factor(NSAIDS, levels = c(0,1), labels = c('No', 'Yes'))
  diabetesT2<-factor(diabetesT2, levels = c(0,1), labels = c('No', 'Yes'))
  hormone_replace<-factor(hormone_replace, levels = c(0,1), labels = c('No', 'Yes'))
  evercol<-factor(evercol, levels = c(0,1), labels = c('No', 'Yes'))
  chemradther<-factor(chemradther, levels = c('Nein', 'Ja'), labels = c('No', 'Yes'))
  neotreat<-factor(neotreat, levels = c(0,1), labels = c('No', 'Yes'))
  
  cardiovad<-factor(ifelse(cardio1_myocar == 1|cardio2_stroke ==1|
                             cardio3_cirheart==1|cardio4_cardinsuffi==1, 'Yes', 
                           ifelse(cardio1_myocar == 0 & cardio2_stroke ==0 &
                                    cardio3_cirheart== 0 & cardio4_cardinsuffi==0, 'No', NA)))
  
  death_crccp<-ifelse(death_all==1&death_crc==1, 1, 
                      ifelse(death_all==1&death_crc==0, 2, 
                             ifelse(is.na(death_crc), NA, 0)))
  
  
  recurr_cp<-ifelse(death_all==1&recurr==1, 1, 
                    ifelse(death_all==1&recurr==0, 2, 
                           ifelse(is.na(recurr), NA, 0)))
  timey<-timeD/365.25
  recurr_timey<-timeD_recurr/365.25
  
  PFS<-ifelse(death_all==1 | recurr==1, 1,
              ifelse(death_all==0&is.na(recurr), NA, 0))
  
  timeD_PFS<-ifelse(recurr==1 & !is.na(recurr), timeD_recurr, 
                    ifelse(death_all==1, timeD,
                           ifelse(death_all==0 & is.na(recurr), NA, timeD)))
  
  timey_PFS<-timeD_PFS/365.25
  
})

summary(blood_clin)


#abstainers or light, moderate, or heavy drinkers, respectively.
# categorize alcohol consumption by sex
w<-subset(blood_clin, Sex == 'Female')
w<-within.data.frame(w, {
  liftimdrink_cat<-factor(ifelse(lifetimealcohol_day ==0, 'Abstainers', 
                                 ifelse(lifetimealcohol_day>0 & lifetimealcohol_day<=12, 'Light drinkers',
                                        ifelse(lifetimealcohol_day>12 & lifetimealcohol_day<=25, 'Moderate drinkers', 
                                               'Heavy drinkers'))),
                          levels = c('Abstainers', 'Light drinkers','Moderate drinkers','Heavy drinkers'), 
                          ordered = T)
  
  recentdrink_cat<-factor(ifelse(lastyralcohol_day ==0, 'Abstainers', 
                                 ifelse(lastyralcohol_day>0 & lastyralcohol_day<=12, 'Light drinkers',
                                        ifelse(lastyralcohol_day>12 & lastyralcohol_day<=25, 'Moderate drinkers', 
                                               'Heavy drinkers'))), 
                          levels = c('Abstainers', 'Light drinkers','Moderate drinkers','Heavy drinkers'), 
                          ordered = T)
})

summary(w$liftimdrink_cat)
summary(w$lastyralcohol_day) # 17 missing
summary(w$recentdrink_cat)

m<-subset(blood_clin, Sex == 'Male')
m<-within.data.frame(m, {
  liftimdrink_cat<-factor(ifelse(lifetimealcohol_day ==0, 'Abstainers', 
                                 ifelse(lifetimealcohol_day>0 & lifetimealcohol_day<=24, 'Light drinkers',
                                        ifelse(lifetimealcohol_day>24 & lifetimealcohol_day<=50, 'Moderate drinkers', 
                                               'Heavy drinkers'))), 
                          levels = c('Abstainers', 'Light drinkers','Moderate drinkers','Heavy drinkers'), 
                          ordered = T)
  
  recentdrink_cat<-factor(ifelse(lastyralcohol_day ==0, 'Abstainers', 
                                 ifelse(lastyralcohol_day>0 & lastyralcohol_day<=24, 'Light drinkers',
                                        ifelse(lastyralcohol_day>24 & lastyralcohol_day<=50, 'Moderate drinkers', 
                                               'Heavy drinkers'))), 
                          levels = c('Abstainers', 'Light drinkers','Moderate drinkers','Heavy drinkers'), 
                          ordered = T)
})

summary(m$liftimdrink_cat)
summary(m$lastyralcohol_day) # 10 missing
summary(m$recentdrink_cat)

blood_clin<-bind_rows(w, m)
summary(blood_clin)

## score calculation
summary(blood_clin$time_blood)
blood_clin<-subset(blood_clin, time_blood>=-10)
blood_clin$time_blood[blood_clin$time_blood<0]<-0 # 2126
save(blood_clin, 
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_clin.RData")

#################################  make table one ###############################
library(tableone)
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_clin.RData")
summary(blood_clin)
dput(names(blood_clin))

vars<-c( "Age_diag", "Sex", "Schooling_years", "Location",
         "TNM_adj", "chemradther",  "neotreat",
         "BMI", "BMI_cat", "BMI_5_14earlier", "BMI5_14er_cat",
         "smoking",  "liftimdrink_cat", # average lifetime daily alcohol consumption
         "physical_activity_lifeav",  
         "cardiovad", "statins", "NSAIDS", "diabetesT2","high_blood_pressure", "hyperlipidemia", "other_cancer")

nonNormalVars<-c("Age_diag", "BMI", "BMI_1_10earlier", "BMI_5_14earlier")

table1<- CreateTableOne(vars = vars,  data = blood_clin, includeNA =T, test = T)
b<-print(table1, nonnormal = nonNormalVars,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$rownames.b.<-substring(total$rownames.b., first = 3)

library(writexl)
write_xlsx(total, path = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/results/Table1_blood.xlsx", col_names = T)

################## score calculation #########
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/beta_ILM_blood.RData")
beta_ILM$uid<-NULL
rownames(beta_ILM)<-beta_ILM$tn_id
beta_ILM$tn_id<-NULL
scores <- read_excel("methy_BMI/BMI_cpgs.xlsx")
View(scores)
View(beta_ILM)
cpg<-subset(scores, CpGs != '(Intercept)')
cpg_orig<-colnames(beta_ILM)
intersect<-Reduce(intersect, list(cpg_orig, cpg$CpGs)) # 4920
scores<-scores[scores$CpGs %in% intersect, ]
summary(as.factor(scores$Study))

##
score<-subset(scores, Study =='Hamilton et al. 2019')
beta_hmn<-beta_ILM[, score$CpGs]

# calculate the score:
for (i in 1:nrow(score)){
  beta_hmn[, i]<-beta_hmn[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_hmn))
beta_hmn$PI_hmn<-rowSums(beta_hmn, na.rm = T)+0.126459884161198
beta_hmn$tn_id<-rownames(beta_hmn)
N<-ncol(beta_hmn)
beta_hmn<-beta_hmn[, c(N, N-1)]
summary(beta_hmn$PI_hmn)

##
score<-subset(scores, Study =='McCartney et al. 2018')
beta_mc<-beta_ILM[, score$CpGs]
for (i in 1:nrow(score)){
  beta_mc[, i]<-beta_mc[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_mc))
beta_mc$PI_mc<-rowSums(beta_mc, na.rm = T)
beta_mc$tn_id<-rownames(beta_mc)
N<-ncol(beta_mc)
beta_mc<-beta_mc[, c(N, N-1)]
summary(beta_mc$PI_mc)

###
score<-subset(scores, Study =='Do et al. 2022')
beta_do<-beta_ILM[, score$CpGs]

# calculate the score:
for (i in 1:nrow(score)){
  beta_do[, i]<-beta_do[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_do))
beta_do$PI_do<-rowSums(beta_do, na.rm = T)+1.45985838513289
beta_do$tn_id<-rownames(beta_do)
N<-ncol(beta_do)
beta_do<-beta_do[, c(N, N-1)]
summary(beta_do$PI_do)

score<-subset(scores, Study =='Mendelson et al. 2017')
beta_mdson<-beta_ILM[, score$CpGs]

# calculate the score:
for (i in 1:nrow(score)){
  beta_mdson[, i]<-beta_mdson[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_mdson))
beta_mdson$PI_mdson<-rowSums(beta_mdson, na.rm = T)
beta_mdson$tn_id<-rownames(beta_mdson)
N<-ncol(beta_mdson)
beta_mdson<-beta_mdson[, c(N, N-1)]
summary(beta_mdson$PI_mdson)

###
score<-subset(scores, Study =='Merzbacher et al. 2023')
beta_mzbach<-beta_ILM[, score$CpGs]

# calculate the score:
for (i in 1:nrow(score)){
  beta_mzbach[, i]<-beta_mzbach[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_mzbach))
beta_mzbach$PI_mzbach<-rowSums(beta_mzbach, na.rm = T)-0.550839855225047
beta_mzbach$tn_id<-rownames(beta_mzbach)
N<-ncol(beta_mzbach)
beta_mzbach<-beta_mzbach[, c(N, N-1)]
summary(beta_mzbach$PI_mzbach)

## merge the three scores together
ilm<-merge(beta_mc, beta_hmn, by = 'tn_id')
ilm<-merge(ilm, beta_do, by = 'tn_id')
ilm<-merge(ilm, beta_mdson, by = 'tn_id')
ilm<-merge(ilm, beta_mzbach, by = 'tn_id')

save(ilm,
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/scores.rdata")

######################## perform multiple impuation ##########################################
############ Blood sample imputation #################
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_clin.RData")

#### start multiple imputation #########
library(mice)
library(survival)
library(mitools)
# calculate cumulative hazard###########change the outcome
dput(names(blood_clin))

clinimput<-subset(blood_clin, select = c("tn_id",  "Age_diag", "Sex", "Schooling_years", "TNM_adj", "chemradther", "neotreat",
                                         "Location", "lifetimealcohol_day", 
                                         "smoking", "physical_activity_lifeav", "hyperlipidemia",
                                         "statins", "NSAIDS", "diabetesT2",  "high_blood_pressure", "cardiovad","other_cancer", 
                                         "height", "BMI", "BMI_5_14earlier", 
                                         "weight_20yr", "weight_30yr", "weight_40yr", "weight_50yr", "weight_60yr", 
                                         "weight_70yr", "weight_80yr",
                                         "timeD", 
                                         "death_all"))
summary(clinimput)

HT1 <- summary(survival::survfit(Surv(timeD, death_all)~1,data=clinimput))
clinimput$haz_os <-  approx(c(0, HT1$time), -log(c(1,HT1$surv)),xout=clinimput$time,method="constant",f=0,rule=2)$y

# set outcome as factor
clinimput$death_all<- factor(clinimput$death_all)

##
# see all the default settings for imputation
impu_default <- mice(clinimput, maxit = 0)
summary(impu_default)

# see the predictor structure
pred <- quickpred(clinimput, exclude = c("id","timey"))
pred

meth <- impu_default$meth
meth

# multiple imputation for 20 times
blood_imputation_20 <- mice(clinimput, maxit = 10, m = 20, seed = 1234, pred = pred, meth = meth, print = TRUE)

save(blood_imputation_20, 
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_imputation_20.RData')

blood_imputatedpro <- vector(20,mode="list")
for (i in 1:20) {
  blood_imputatedpro[[i]] <- mice::complete(blood_imputation_20, i)
  blood_imputatedpro[[i]]$haz_os<-NULL
  blood_imputatedpro[[i]]$death_all<-NULL
  blood_imputatedpro[[i]]$timeD<-NULL
}

summary(blood_imputatedpro[[1]])
save(blood_imputatedpro,
     file = 'C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_imputatedpro.RData')

###### calculate the WYOs and cumulative BMI in the 20 imputed dataset ####
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_imputatedpro.RData")
dput(names(df))

BMIcum<-vector(20, mode = "list")
for (i in 1:20) {
  blood_imputatedpro[[i]]$height<-blood_imputatedpro[[i]]$height/100
  BMIcum[[i]] <- subset(blood_imputatedpro[[i]],
                        select = c("tn_id", "Age_diag",   "height", "BMI", "weight_20yr", 
                                   "weight_30yr", "weight_40yr", "weight_50yr", "weight_60yr", "weight_70yr", 
                                   "weight_80yr"))
  
  BMIcum[[i]]$Year_recent<-BMIcum[[i]]$Age_diag -5 ## account for weigthloss due to tumor
  
  BMIcum[[i]]$BMI_20 <- BMIcum[[i]]$weight_20yr/(BMIcum[[i]]$height^2)
  BMIcum[[i]]$BMI_30 <- BMIcum[[i]]$weight_30yr/(BMIcum[[i]]$height^2)
  BMIcum[[i]]$BMI_40 <- BMIcum[[i]]$weight_40yr/(BMIcum[[i]]$height^2)
  BMIcum[[i]]$BMI_50 <- BMIcum[[i]]$weight_50yr/(BMIcum[[i]]$height^2)
  BMIcum[[i]]$BMI_60 <- BMIcum[[i]]$weight_60yr/(BMIcum[[i]]$height^2)
  BMIcum[[i]]$BMI_70 <- BMIcum[[i]]$weight_70yr/(BMIcum[[i]]$height^2)
  BMIcum[[i]]$BMI_80 <- BMIcum[[i]]$weight_70yr/(BMIcum[[i]]$height^2)
  
  BMIcum[[i]] <- subset(BMIcum[[i]],
                        select = c("tn_id", "Age_diag", "Year_recent",  "BMI",  
                                   "BMI_20", "BMI_30", "BMI_40", "BMI_50", "BMI_60", 
                                   "BMI_70", "BMI_80"))

}

summary(BMIcum[[1]])
save(BMIcum, 
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_BMIcum20raw.RData")

load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_BMIcum20raw.RData")

# Initialize an empty list to store results
cumBMI_result_list <- list()

# Loop through each dataset in BMIcum
for (dataset_index in 1:20) {
  # Load the current dataset (e.g., BMIcum[[1]], BMIcum[[2]], ..., BMIcum[[20]])
  df <- BMIcum[[dataset_index]]
  
  # Convert from wide to long format
  long_data <- reshape(
    df,
    varying = list(c("BMI_20", "BMI_30", "BMI_40", "BMI_50", "BMI_60", "BMI_70", "BMI_80")), 
    v.names = "BMI_dyna",
    timevar = "Age_measure",
    times = c(20, 30, 40, 50, 60, 70, 80),
    direction = "long"
  )
  
  long_data <- long_data[order(long_data$tn_id), ]
  long_data <- long_data %>% 
    group_by(tn_id) %>%
    filter(Age_measure <= Age_diag)
  
  # Initialize an empty list to store results for each dataset
  result_list <- list()
  
  # Get unique 'id' values
  unique_ids <- unique(long_data$id)
  
  # Iterate over each 'id'
  for (i in unique_ids) {
    subset_df <- subset(long_data, id == i)
    
    # Perform linear interpolation
    interpolation <- approx(
      x = subset_df$Age_measure,
      y = subset_df$BMI_dyna,
      xout = seq(20, unique(subset_df$Age_diag)),
      rule = 2
    )
    
    group_result <- data.frame(
      tn_id = unique(subset_df$tn_id),
      Age_diag = unique(subset_df$Age_diag),
      Year_recent = unique(subset_df$Year_recent),
      Age_measure = interpolation$x,
      BMI_dyna = interpolation$y
    )
    
    # Append the result to the list for this dataset
    result_list[[i]] <- group_result
  }
  
  # Combine all results for this dataset into a single dataframe
  final_result <- do.call(rbind, result_list)
  final_result$eBMI_dyn <- final_result$BMI_dyna - 25
  final_result$eBMI_dyn[final_result$eBMI_dyn < 0] <- 0
  
  # Calculate WYO, cuBMI, etc.
  
  # 1. WYO_tdiagnosis: eBMI from 20 to diagnosis
  WYO_td <- aggregate(eBMI_dyn ~ tn_id, 
                      data = final_result, FUN = sum, na.rm = TRUE)
  colnames(WYO_td)[2] <- 'WYO_td'
  
  # 2. WYO_tr5yr: eBMI from 20 to 5-14 years before diagnosis
  final_result$lower_tens <- floor(final_result$Year_recent / 10) * 10
  BMI_514 <- final_result %>% 
    group_by(tn_id) %>%
    filter(Age_measure <= lower_tens)
  
  WYO_tr5yr <- aggregate(eBMI_dyn ~ tn_id, 
                         data = BMI_514, FUN = sum, na.rm = TRUE)
  colnames(WYO_tr5yr)[2] <- 'WYO_tr5yr'
  
  # 3. Median and Mean BMI from 20 to diagnosis / 5-14 years before diagnosis
  cuBMI_mediandiagn <- aggregate(BMI_dyna ~ tn_id, 
                                 data = final_result, FUN = median, na.rm = TRUE)
  
  cuBMI_meandiagn <- aggregate(BMI_dyna ~ tn_id, 
                               data = final_result, FUN = mean, na.rm = TRUE)
  
  cuBMI_medianearl <- aggregate(BMI_dyna ~ tn_id, 
                                data = BMI_514, FUN = median, na.rm = TRUE)
  
  cuBMI_meanearl <- aggregate(BMI_dyna ~ tn_id, 
                              data = BMI_514, FUN = mean, na.rm = TRUE)
  
  colnames(cuBMI_meandiagn)[2] <- 'meanBMI_tdiag'
  colnames(cuBMI_mediandiagn)[2] <- 'medianBMI_tdiag'
  colnames(cuBMI_medianearl)[2] <- 'medianBMI_t514'
  colnames(cuBMI_meanearl)[2] <- 'meanBMI_t514'
  
  # Bind all results together for this dataset
  cumBMI_result <- merge(WYO_td, WYO_tr5yr, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_mediandiagn, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_meandiagn, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_medianearl, by = 'tn_id')
  cumBMI_result <- merge(cumBMI_result, cuBMI_meanearl, by = 'tn_id')
  
  # Store the result of this dataset in the list
  cumBMI_result_list[[dataset_index]] <- cumBMI_result
}

d<-cumBMI_result_list[[1]]
summary(d)
save(cumBMI_result_list, 
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/cumBMI_result_list20.RData")

########### merge cumBMI and baseline variablesa and methylation-based BMI #####
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/cumBMI_result_list20.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_imputatedpro.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_clin.RData")
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/scores.rdata")
dput(names(blood_clin))
df_clin<-subset(blood_clin, 
                select = c("tn_id", "lab_nummer", 
                           "late_days",  "time_blood",
                           "timeD", "timeM", "death_all", "death_crc", "timeD_recurr", "timeM_recurr", 
                           "recurr", "recurr_type", "recurr_cp", "death_crccp"))

colnames(ilm)[1]<-'lab_nummer'

df_clin<-merge(df_clin, ilm, by = 'lab_nummer', all.x = T)

blood_all<-vector(20, mode = "list")

for (i in 1:20) {
  blood_imputatedpro[[i]]<-subset(blood_imputatedpro[[i]],
                                  select =  c("tn_id","Age_diag", "Sex", "Schooling_years", "TNM_adj", "chemradther", 
                                              "neotreat", "Location", "lifetimealcohol_day", "smoking", "physical_activity_lifeav", 
                                              "hyperlipidemia", "statins", "NSAIDS", "diabetesT2", "high_blood_pressure", 
                                              "cardiovad", "other_cancer", "BMI", "BMI_5_14earlier", "height"))

  
  blood_all[[i]]<-merge(blood_imputatedpro[[i]],
                        cumBMI_result_list[[i]], by = 'tn_id')
  
  blood_all[[i]]<-merge(blood_all[[i]], df_clin, by = 'tn_id')
  
}

summary(blood_all[[1]])
save(blood_all,
     file = "C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")



