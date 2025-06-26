### use another preproccessing strategy
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/bloodm_ILM.RData")
scores <- read_excel("alchohol_methy/scores.xlsx", 
                     sheet = "alcohol_socre_details")
cpg<-as.data.frame(rownames(dataset_ILM))


intersect<-Reduce(intersect, list(scores$CpG, cpg$`rownames(dataset_ILM)`)) # 588
beta_ILM<-dataset_ILM[intersect, ]

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

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/blood_methylation/processed_data/methy_id.rdata")

beta_ILM<-merge(id, beta_ILM, by = 'uid')
beta_ILM[intersect]<-lapply(beta_ILM[intersect], function(x)as.numeric(x))
save(beta_ILM,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/beta_ILM.RData")

### score calculation
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/beta_ILM.RData")
beta_ILM$uid<-NULL
rownames(beta_ILM)<-beta_ILM$tn_id
beta_ILM$tn_id<-NULL
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
beta_mc$lab_nummer<-rownames(beta_mc)
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
beta_epitob$lab_nummer<-rownames(beta_epitob)
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
beta_liu$lab_nummer<-rownames(beta_liu)
N<-ncol(beta_liu)
beta_liu<-beta_liu[, c(N, N-1)]
summary(beta_liu$PI_liu)

## same results 
## merge the three scores together
ilm<-merge(beta_mc, beta_epitob, by = 'lab_nummer')
ilm<-merge(ilm, beta_liu, by = 'lab_nummer')
## merge with clinical variables
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_clin.RData")
summary(blood_clin)
summary(ilm)
blood_clin<-blood_clin[, c(1:46)]
blood_clin<-merge(blood_clin, ilm, by = 'lab_nummer')
summary(blood_clin)
save(blood_clin,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/dcheck_blood_clin.RData")

#################### plot #######################
load( "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/dcheck_blood_clin.RData")









  
  
