####
library(rms)
library(ggplot2)
library(plyr)
library(ggpubr)
load("C:/Users/T532N/Desktop/BMI_m/methy_BMI/processed_data/blood_allv20.RData")
df<-blood_all[[1]]
summary(df)
## randomly selected one dataset to plot the curves
dd <- datadist(df)
options(datadist="dd")

AIC <- Inf  # Initialize AIC with a high value
nk <- 3     # Default nk value

for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all) ~ rcs(BMI, i),
             data = df, x = TRUE)
  tmp <- AIC(fit)  # Use AIC function from rms package
  if (tmp < AIC) {
    AIC <- tmp
    nk <- i
  }
}

print(nk)
dput(names(df))
### plot the association with CSC
df$TNM_adj <- as.factor(as.character(df$TNM_adj))
is.ordered(df$TNM_adj)
cont_var <- c("BMI", "PI_mdson",  "PI_mc", "PI_hmn", "PI_do",  "PI_mzbach")
adjust_covariates <- "+Age_diag + Sex + Schooling_years + TNM_adj + Location + chemradther +
  neotreat+ smoking+physical_activity_lifeav + lifetimealcohol_day + cardiovad + NSAIDS +
  statins + diabetesT2 + high_blood_pressure + hyperlipidemia+other_cancer"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each BMI measure
for (var in cont_var) {
  
  # Fit the model
  formula <- as.formula(paste("Surv(time = time_blood, time2 = timeD, death_all) ~ rcs(", var, ", 3)", adjust_covariates))
  fit <- cph(formula, data = df, x = TRUE)
  
  # Generate predictions
  surv2 <- Predict(fit, name = var, ref.zero = TRUE, fun = exp)
  
  # Calculate quartiles for the BMI variable
  t <- quantile(df[[var]], probs = seq(0, 1, 1/4))[c(2, 3, 4)]
  
  # Create the ggplot for the current BMI measure
  plot <- ggplot(surv2,adj.subtitle = FALSE, aes(x = !!sym(var), y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "pink") +
    geom_line(color = "red") +
    labs(x = var, y = "Hazard ratio", title = paste(var, "(OS)")) +
    geom_vline(xintercept = t, linetype = "dashed", color = "gray") +
    coord_cartesian(ylim = c(0, 2.5)) +  # Use this instead of ylim()
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "top",
          legend.spacing.x = unit(1.2, 'cm'),
          legend.text = element_text(size = 11),
          axis.line.x = element_line(color = "black", size = 0.8),
          axis.line.y = element_line(color = "black", size = 0.8),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm"))) +
    annotate("text", x = t[1], y = 0, label = "Q1", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    annotate("text", x = t[2], y = 0, label = "Q2", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    annotate("text", x = t[3], y = 0, label = "Q3", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic")
  
  # Add the plot to the list
  plot_list[[var]] <- plot
}

plot_list[[3]]

os<-ggarrange(plot_list[[1]], plot_list[[2]],plot_list[[3]],
          plot_list[[4]], plot_list[[5]],plot_list[[6]], nrow = 2, ncol = 3)

## 6 *11 

##### cancer-specific survival ###
df$TNM_adj <- as.factor(as.character(df$TNM_adj))
is.ordered(df$TNM_adj)
cont_var <- c("BMI", "PI_mdson",  "PI_mc", "PI_hmn", "PI_do",  "PI_mzbach")
adjust_covariates <- "+Age_diag + Sex + Schooling_years + TNM_adj + Location + chemradther +
  neotreat+ smoking+physical_activity_lifeav + lifetimealcohol_day + cardiovad + NSAIDS +
  statins + diabetesT2 + high_blood_pressure + hyperlipidemia+other_cancer"

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each BMI measure
for (var in cont_var) {
  
  # Fit the model
  formula <- as.formula(paste("Surv(time = time_blood, time2 = timeD, death_crccp==1) ~ rcs(", var, ", 3)", adjust_covariates))
  fit <- cph(formula, data = df, x = TRUE)
  
  # Generate predictions
  surv2 <- Predict(fit, name = var, ref.zero = TRUE, fun = exp)
  
  # Calculate quartiles for the BMI variable
  t <- quantile(df[[var]], probs = seq(0, 1, 1/4))[c(2, 3, 4)]
  
  # Create the ggplot for the current BMI measure
  plot <- ggplot(surv2,adj.subtitle = FALSE, aes(x = !!sym(var), y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "pink") +
    geom_line(color = "red") +
    labs(x = var, y = "Hazard ratio", title = paste(var, "(CSS)")) +
    geom_vline(xintercept = t, linetype = "dashed", color = "gray") +
    coord_cartesian(ylim = c(0, 2.5)) +  # Use this instead of ylim()
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "top",
          legend.spacing.x = unit(1.2, 'cm'),
          legend.text = element_text(size = 11),
          axis.line.x = element_line(color = "black", size = 0.8),
          axis.line.y = element_line(color = "black", size = 0.8),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm"))) +
    annotate("text", x = t[1], y = 0, label = "Q1", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    annotate("text", x = t[2], y = 0, label = "Q2", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic") +
    annotate("text", x = t[3], y = 0, label = "Q3", color = '#8d99ae', family = "Arial", size = 4, fontface = "italic")
  
  # Add the plot to the list
  plot_list[[var]] <- plot
}

plot_list[[3]]

css<-ggarrange(plot_list[[1]], plot_list[[2]],plot_list[[3]],
               plot_list[[4]], plot_list[[5]],plot_list[[6]], nrow = 2, ncol = 3)

drplot<-ggarrange(os, css, nrow = 2)









