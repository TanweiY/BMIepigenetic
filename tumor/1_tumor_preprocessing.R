
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

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/methy_normILM.RData")
beta_ILM<-impute.knn(dataset_ILM)$data

save(beta_ILM, 
     file = "/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/tumorbeta_ILM_imputed.RData")

scores <-  read_excel("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/methy_BMI/BMI_cpgs.xlsx")

cpg<-as.data.frame(rownames(beta_ILM)) 

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
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/methy_id.RData")

dataset_ILM<-merge(id, dataset_ILM, by = 'uid')
dataset_ILM[intersect]<-lapply(dataset_ILM[intersect], function(x)as.numeric(x))

save(beta_ILM,
     file = "/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/beta_ILM_tumor.RData")


## merge with clinical variables
library(haven)
baseline<-read_sas('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/dachs_cohort/dachskomplett_categ_20191231.sas7bdat')
variables <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/variables.xlsx", 
                        sheet = "baseline")

base_var<-c(variables$var)
base<-subset(baseline, select = base_var)
colnames(base)<-c(variables$Names)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/dataset_ILM_tumor.RData")

base<-base[(base$tn_id %in% dataset_ILM$id), ] ## 2279   
rm(dataset_ILM)

## select follow up dataset
follow<-read_sas('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/dachs_cohort/follow.sas7bdat')
variables <- read_excel("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/variables.xlsx", 
                        sheet = "follow")
f<-subset(follow, select = variables$var)
colnames(f)<-variables$Names

tumor_clin<-merge(base, f, by = 'tn_id')

## exclude without lifetime alchohol consumption
tumor_clin<-subset(tumor_clin, !is.na(tumor_clin$lifetimealcohol_day))

## variable coding ##
## N = 2279
summary(tumor_clin)
colnames(tumor_clin)
tumor_clin<-within.data.frame(tumor_clin, {
  time_blood<-as.numeric(difftime(Blutabnahme_Dat, indexdat,  units = "days"))
  BMI_cat<- cut(BMI, breaks = c(-Inf, 18.5, 25, 30, Inf), labels = c('<18.5','18.5-25','25-30', '≥30'))
  BMIearlier_cat<- cut(BMI_earlier, breaks = c(-Inf, 18.5, 25, 30, Inf), labels = c('<18.5','18.5-25','25-30', '≥30'))
  Schooling_years<- factor(Schooling_years, levels = c(1,2,3), labels = c('<9','9-10','>10'))
  
  Sex<-factor(Sex, levels = c(1,2), labels = c('Female', 'Male'))
  TNM_adj<-factor(TNM_adj, levels = c(1,2,3,4), labels = c('I', 'II', 'III', 'IV'), ordered = T)
  smoking<-factor(smoking, levels = c(0, 1, 2), labels = c('Never', 'Former', 'Current'))
  statins<-factor(statins, levels = c(0,1), labels = c('No', 'Yes'))
  other_cancer<-factor(other_cancer, levels = c(0,1), labels = c('No', 'Yes'))
  high_blood_pressure<-factor(high_blood_pressure, levels = c(0,1), labels = c('No', 'Yes'))
  NSAIDS<-factor(NSAIDS, levels = c(0,1), labels = c('No', 'Yes'))
  diabetesT2<-factor(diabetesT2, levels = c(0,1), labels = c('No', 'Yes'))
  hormone_replace<-factor(hormone_replace, levels = c(0,1), labels = c('No', 'Yes'))
  evercol<-factor(evercol, levels = c(0,1), labels = c('No', 'Yes'))
  chemradther<-factor(chemradther, levels = c('Nein', 'Ja'), labels = c('No', 'Yes'))
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
})

summary(tumor_clin)

### classification of alchohol consumption
# Women with consumption levels of 0, .0– 12, .12–25, or .25 g/d and
#men with consumption levels of 0, .0–24, .24–50, or .50 g/d were each classified as 
w<-subset(tumor_clin, Sex == 'Female')
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

m<-subset(tumor_clin, Sex == 'Male')
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


tumor_clin<-bind_rows(w, m)
summary(tumor_clin)
save(tumor_clin,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")

################ make table one #############################
library(tableone)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")

vars<-c( "Age_diag", "Sex", "Schooling_years","liftimdrink_cat","recentdrink_cat",
         "TNM_adj", "chemradther", "crc2sites","BMI_cat",  "BMIearlier_cat",
         "smoking", "physical_activity_lifeav", "evercol", "hormone_replace",
         "cardiovad", "statins", "NSAIDS", "diabetesT2","high_blood_pressure","other_cancer")

nonNormalVars<-c("Age_diag")
table1<- CreateTableOne(vars = vars,  data = tumor_clin, includeNA =T, test = T)

b<-print(table1, nonnormal = nonNormalVars,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$rownames.b.<-substring(total$rownames.b., first = 3)

library(writexl)
write_xlsx(total, path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/Table1_tumor.xlsx", col_names = T)

#################### score calculation ###################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/dataset_ILM_tumor.RData")
beta_ILM<-dataset_ILM
beta_ILM$uid<-NULL
rownames(beta_ILM)<-beta_ILM$id
beta_ILM$id<-NULL
library(readxl)
scores <- read_excel("alchohol_methy/scores.xlsx", sheet = "alcohol_socre_details")
View(scores)
View(beta_ILM)
cpg<-scores$CpG[-1]
cpg_orig<-colnames(beta_ILM)
intersect<-Reduce(intersect, list(cpg_orig, cpg))
scores<-scores[scores$CpG %in% intersect, ]
score<-subset(scores, Study =='McCartney et al. 2018')
beta_mc<-beta_ILM[, score$CpG]
# calculate the score:
for (i in 1:nrow(score)){
  beta_mc[, i]<-beta_mc[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_mc))
beta_mc$PI_mc<-rowSums(beta_mc, na.rm = T)
View(beta_mc)
beta_mc$tn_id<-rownames(beta_mc)
N<-ncol(beta_mc)
beta_mc<-beta_mc[, c(N, N-1)]
summary(beta_mc$PI_mc)

score<-subset(scores, Study =='Chamberlain et al. 2022')
beta_epitob<-beta_ILM[, score$CpG]
for (i in 1:nrow(score)){
  beta_epitob[, i]<-beta_epitob[, i]*as.numeric(score[i, 3])
}
sum(is.na(beta_epitob))
beta_epitob$PI_epitob<-rowSums(beta_epitob, na.rm = T)+94.583
beta_epitob$tn_id<-rownames(beta_epitob)
N<-ncol(beta_epitob)
beta_epitob<-beta_epitob[, c(N, N-1)]
summary(beta_epitob$PI_epitob)

## the 3rd score ##
library(dnamalci)
score<-subset(scores, Study =='Liu et al. 2018')
beta_liu<-beta_ILM[, score$CpG]
matrix_liu<-t(beta_liu)
dnam.alc <- dnamalci(matrix_liu)
colnames(matrix_liu)

## double check with another method to calculate
coef<-as.data.frame(dnam.alc$coefficients)
coef$CpG<-rownames(coef)
colnames(coef)[1]<-'coefficients'

beta_liu<-beta_ILM[, coef$CpG]
# calculate the coef:
for (i in 1:nrow(coef)){
  beta_liu[, i]<-beta_liu[, i]*as.numeric(coef[i, 1])
}
sum(is.na(beta_liu))
beta_liu$PI_liu<-rowSums(beta_liu, na.rm = T)
View(beta_liu)
beta_liu$tn_id<-rownames(beta_liu)
N<-ncol(beta_liu)
beta_liu<-beta_liu[, c(N, N-1)]
summary(beta_liu$PI_liu)

## merge the three scores together
ilm<-merge(beta_mc, beta_epitob, by = 'tn_id')
ilm<-merge(ilm, beta_liu, by = 'tn_id')

## merge with clinical variables
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")
tumor_clin<-merge(tumor_clin, ilm, by = 'tn_id')
summary(tumor_clin)
save(tumor_clin,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")


















