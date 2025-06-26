# 5.17 5 year follow up
#"How many day per week did you usually drink alcohol in the past 12 months?" 
# ls02alktage 
# ls02abier bottles/week
# ls02awein glass/week
# ls02aschnaps glass/week
# Does the participant currently drink alcohol at least one day per week?  _alktrinker

# 6.8 10 year follow-up
# z_ls02alktage  # On how many days per week have you drank alcohol in the past 12 months (no alcohol=0 days)?
# z_ls02abier
# z_ls02awein
# z_ls02aschnaps

follow<-read_sas('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/dachs_cohort/follow.sas7bdat')
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/blood_clin.RData")
f_dr<-follow[(follow$tn_id %in% blood_clin$tn_id), ]
rm(follow)
var<-c('tn_id', 'ls02alktage', 'ls02abier', 'ls02awein', 'ls02aschnaps', 
       'z_ls02alktage', 'z_ls02bier', 'z_ls02wein', 'z_ls02schnaps')

f_dr<-subset(f_dr, select = var)

sex<-blood_clin[, c(2, 7)]
f_dr<-merge(sex, f_dr, by = 'tn_id')

summary(f_dr)

## make the variable myself
# Constants for ethanol content in g/100 mL
ethanol_beer <- 4    # g/100mL
ethanol_wine <- 8.6  # g/100mL
ethanol_liquor <- 33 # g/100mL

# Serving sizes
serving_beer_ml <- 330 # ml per bottle
serving_wine_ml <- 250 # ml per glass
serving_liquor_ml <- 20 # ml per glass

# Function to calculate ethanol consumption
calculate_ethanol <- function(beer_week, wine_week, liquor_week) {
  # Weekly ethanol consumption in grams
  beer_ethanol_week <- (beer_week * serving_beer_ml * ethanol_beer) / 100
  wine_ethanol_week <- (wine_week * serving_wine_ml * ethanol_wine) / 100
  liquor_ethanol_week <- (liquor_week * serving_liquor_ml * ethanol_liquor) / 100
  
  # Total daily ethanol consumption in grams
  total_ethanol_daily <- (beer_ethanol_week + wine_ethanol_week + liquor_ethanol_week) / 7
  
  return(total_ethanol_daily)
}

## function to categorize drinkers #####
categorize_drinker <- function(ethanol_daily, gender) {
  # Check if either ethanol_daily or gender is NA
  if (is.na(ethanol_daily) || is.na(gender)) {
    return(NA)
  }
  
  # Classify based on gender and ethanol consumption
  if (gender == "Female") {
    if (ethanol_daily == 0) {
      return("Abstainer")
    } else if (ethanol_daily <= 12) {
      return("Light")
    } else if (ethanol_daily <= 25) {
      return("Moderate")
    } else {
      return("Heavy")
    }
  } else if (gender == "Male") {
    if (ethanol_daily == 0) {
      return("Abstainer")
    } else if (ethanol_daily <= 24) {
      return("Light")
    } else if (ethanol_daily <= 50) {
      return("Moderate")
    } else {
      return("Heavy")
    }
  }
}

# Apply the function to calculate total daily ethanol consumption
f_dr$total5_ethanol_daily <- mapply(calculate_ethanol, f_dr$ls02abier, f_dr$ls02awein, f_dr$ls02aschnaps)
summary(f_dr$total5_ethanol_daily)
summary(f_dr$ls02alktage)
f_dr$total5_ethanol_daily[f_dr$ls02alktage == 0]<-0
f_dr$drinking_5yr <- mapply(categorize_drinker, f_dr$total5_ethanol_daily, f_dr$Sex)
summary(as.factor(f_dr$drinking_5yr))

f_dr$total10_ethanol_daily <- mapply(calculate_ethanol, f_dr$z_ls02bier, f_dr$z_ls02wein, f_dr$z_ls02schnaps)
f_dr$total10_ethanol_daily[f_dr$z_ls02alktage == 0]<-0
f_dr$drinking_10yr <- mapply(categorize_drinker, f_dr$total10_ethanol_daily, f_dr$Sex)
summary(as.factor(f_dr$drinking_10yr))

save(f_dr, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/f_drink.RData")

## start table one ##
library(tableone)

dput(names(f_dr))

vars<-c("drinking_5yr", "drinking_10yr")

table1<- CreateTableOne(vars = vars,  data = f_dr, includeNA =F, test = T)

b<-print(table1,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$rownames.b.<-substring(total$rownames.b., first = 3)

library(writexl)
write_xlsx(total, path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/results/drink_follow.xlsx", col_names = T)





