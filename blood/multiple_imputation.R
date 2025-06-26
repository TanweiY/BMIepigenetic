## 
############ Blood sample imputation #################
load("/omics/odcf/analysis/OE0167_projects_temp/dachs_genetic_data_platform/methy_BMI/processed_data/blood_clin.RData")

#### start multiple imputation #########
library(mice)
library(survival)
library(mitools)
# calculate cumulative hazard###########change the outcome
summary(blood_clin$hormone_replace)

clinimput<-subset(blood_clin, select = c("tn_id",  "indexjahr", 
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

save(clinimput, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/bloodclinimpute_orig.RData')

load('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/bloodclinimpute_orig.RData')
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

save(blood_imputation_20, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_imputation_20.RData')

####
######### preprocess the imputed results #######
load('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_imputation_20.RData')
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_clin.RData")
## subset the results based on sample source
dput(names(blood_clin))

condata<-subset(blood_clin, 
                select = c("tn_id", "time_blood", "timeD","death_all", "death_crccp","timeD_recurr","recurr_cp",
                           "recentdrink_cat",  "liftimdrink_cat", "PI_mc", "PI_epitob", "PI_liu", "PI_Liu_md", 
                           "PI_mc_md", "PI_epitob_md"))

blood_imputatedpro <- vector(20,mode="list")

score<-c('PI_epitob', 'PI_mc', 'PI_liu')

for (i in 1:20) {
  blood_imputatedpro[[i]] <- mice::complete(blood_imputation_20, i)
  blood_imputatedpro[[i]]$haz_os<-NULL
  blood_imputatedpro[[i]]$death_all<-NULL
  blood_imputatedpro[[i]]$timeD<-NULL
  blood_imputatedpro[[i]]<-merge(blood_imputatedpro[[i]], condata, by = 'tn_id')
  blood_imputatedpro[[i]][score]<-scale(blood_imputatedpro[[i]][score])
  
}



summary(blood_imputatedpro[[1]])
save(blood_imputatedpro,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_imputatedpro.RData')











