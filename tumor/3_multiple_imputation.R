## 
############ Tumor sample imputation #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.RData")

#### start multiple imputation #########
library(mice)
library(survival)
library(mitools)
# calculate cumulative hazard###########change the outcome
summary(tumor_clin$hormone_replace)

clinimput<-subset(tumor_clin, select = c("tn_id",  "indexjahr", 
                                         "Age_diag", "Sex", "Schooling_years", "TNM_adj", "crc2sites", 
                                         "lastyralcohol_day", "lifetimealcohol_day", "BMI", "BMI_earlier", 
                                         "smoking", "physical_activity_lifeav", "evercol", "hormone_replace", 
                                         "statins", "NSAIDS", "diabetesT2",  "high_blood_pressure", "cardiovad",
                                         "other_cancer", "chemradther", "neotreat",  "timeD", 
                                         "death_all"))
summary(clinimput)

HT1 <- summary(survival::survfit(Surv(timeD, death_all)~1,data=clinimput))
clinimput$haz_os <-  approx(c(0, HT1$time), -log(c(1,HT1$surv)),xout=clinimput$time,method="constant",f=0,rule=2)$y

# set outcome as factor
clinimput$death_all<- factor(clinimput$death_all)

##
(sum(is.na(clinimput))/prod(dim(clinimput)))*100 --> 0.4% missing. 

# see all the default settings for imputation
impu_default <- mice(clinimput, maxit = 0)
summary(impu_default)

# see the predictor structure
pred <- quickpred(clinimput, exclude = c("id","timey"))
pred

meth <- impu_default$meth
meth

# multiple imputation for 20 times
tumor_imputation_20 <- mice(clinimput, maxit = 10, m = 20, seed = 1234, pred = pred, meth = meth, print = TRUE)

save(tumor_imputation_20, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_imputation_20.RData')

####
######### preprocess the imputed results #######
## subset the results based on sample source
load('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_imputation_20.RData')
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")
summary(tumor_clin$liftimdrink_cat)
tumor_clin$liftimdrink_cat<-factor(as.character(tumor_clin$liftimdrink_cat), levels = c('Abstainers', 'Light drinkers', 
                                                                                        'Moderate drinkers',
                                                                                        'Heavy drinkers'))

tumor_clin$recentdrink_cat<-factor(as.character(tumor_clin$recentdrink_cat), levels = c('Abstainers', 'Light drinkers', 
                                                                                        'Moderate drinkers',
                                                                                        'Heavy drinkers'))

save(tumor_clin,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")

condata<-subset(tumor_clin, 
                select = c("tn_id", "late_days", "timeD","death_all", "death_crccp","timeD_recurr","recurr_cp",
                           "recentdrink_cat",  "liftimdrink_cat", "PI_mc", "PI_epitob", "PI_liu"))

tumor_imputatedpro <- vector(20,mode="list")

score<-c('PI_epitob', 'PI_mc', 'PI_liu')

for (i in 1:20) {
  tumor_imputatedpro[[i]] <- mice::complete(tumor_imputation_20, i)
  tumor_imputatedpro[[i]]$haz_os<-NULL
  tumor_imputatedpro[[i]]$death_all<-NULL
  tumor_imputatedpro[[i]]$timeD<-NULL
  tumor_imputatedpro[[i]]<-merge(tumor_imputatedpro[[i]], condata, by = 'tn_id')
  tumor_imputatedpro[[i]][score]<-scale(tumor_imputatedpro[[i]][score])
  
  
}

summary(tumor_imputatedpro[[1]])
save(tumor_imputatedpro,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_imputatedpro.RData')

load()






