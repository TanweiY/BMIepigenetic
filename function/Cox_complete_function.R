### survival function
multicoxfull<-function(data, entrytime, score, outcome){
       if (!outcome %in% c('OS', 'DOC', 'CSS', 'TTR')){
         stop("Invalid outcome. Please choose 'OS', 'DOC', 'TTR', or 'CSS'.")
       }
  
       fml_base <- paste0("Surv(time = ", entrytime, 
                         if (outcome == 'OS') ", time2 = timeD, death_all)" else
                           if (outcome =='DOC') ", time2 = timeD, death_crccp==2)" else 
                             if (outcome =='TTR') ", time2 = timeD_recurr, recurr_cp==1)" else 
                               ", time2 = timeD, death_crccp==1)", 
                         "~Age_diag + Sex + Schooling_years + TNM_adj + Location + chemradther +
                    neotreat+ smoking+physical_activity_lifeav + lifetimealcohol_day + cardiovad + NSAIDS +
                    statins + diabetesT2 + high_blood_pressure + hyperlipidemia+other_cancer+ ")       
       
       fml <- as.formula(paste0(fml_base, score))               
         
       model <- coxph(fml,data = data)
         
       b<-as.data.frame(ShowRegTable(model, exp = TRUE, digits = 2, pDigits = 3,
                                     printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T)
       
       
       score_rows <- grep(score, rownames(b))
       
       b<-b[score_rows, ]
       colnames(b)[c(1, 2)]<-c('aHR(95% CI)','aPvalue')
       b$score<-rownames(b)
       b$outcome<-outcome
       return(b)
}






