##### linear regression
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_imputatedpro.RData")
df<-blood_imputatedpro[[1]]
summary(df)
df<-within.data.frame(df,{
  lifedrinbinary<- factor(ifelse(liftimdrink_cat == 'Abstainers', 'non-drinker','drinker'), 
                          levels = c('non-drinker', 'drinker'))
  
  recentdrinbinary<-  factor(ifelse(recentdrink_cat == 'Abstainers', 'non-drinker',
                                    ifelse(is.na(recentdrink_cat), NA, 'drinker')),
                             levels = c('non-drinker', 'drinker'))
})

dput(names(df))

fml<-formula(paste("lifedrinbinary",  "~", "PI_liu", "+Age_diag + Sex + BMI_earlier + smoking + NSAIDS + hormone_replace + evercol+ physical_activity_lifeav"))
model <- glm(fml, data = df, family = "binomial")
result <- as.data.frame(coef(summary(model)))[2, ]
result <-within.data.frame(result, {
  OR<-round(exp(Estimate), 2)
  LL<-round(exp(confint.default(model))[2], 2)
  UL<-round(exp(confint.default(model))[12], 2)
  #UL<-round(exp(confint.default(model))[6], 2)
  aORCI<-paste0(OR, ' (', LL, ', ', UL, ')')
})


multilog<-function(data, score, outcome){
  fml<-formula(paste(outcome, "~", score, "+Age_diag + Sex + recbmi + alcoholdayweek + nsaidsreg + hrt + evercol+metav"))
  #fml<-formula(paste(outcome, "~", score, "+Age_diag + Sex"))
  model <- glm(fml, data = data, family = "binomial")
  
  result <- as.data.frame(coef(summary(model)))[2, ]
  result <-within.data.frame(result, {
    OR<-round(exp(Estimate), 2)
    LL<-round(exp(confint.default(model))[2], 2)
    UL<-round(exp(confint.default(model))[12], 2)
    #UL<-round(exp(confint.default(model))[6], 2)
    aORCI<-paste0(OR, ' (', LL, ', ', UL, ')')
  })
  
  result<-result[, c(5, 4)]
  colnames(result)[c(1, 2)]<-c('aOR(95% CI)','aPvalue')
  result[, 2]<-round(result[, 2], 3)
  result[, 2]<-ifelse(result[, 2] ==0, '<0.0001', result[, 2])
  result$marker<- rownames(result)
  result$comp<- outcome
  output<-result
  return(output)
}