library(pROC)
get_aucci<-function(dataset, outcome_col, score_col) {
  outcome<-dataset[[outcome_col]]
  score<-dataset[[score_col]]
  r<-roc(outcome, score)
  ci<-ci.auc(r, method = c('bootstrap'), conf.level = 0.95)
  result<-data.frame(Score = score_col, 
                     Outcome = outcome_col,
                     AUCCI = paste0(round(r$auc, 2), ' (', round(ci[1], 2), ', ',  round(ci[2], 2), ')'))
  return (result)
  
}

