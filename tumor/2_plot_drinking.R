library(ggplot2)
library(ggpubr)
library(cgwtools)
library(car)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/alchohol_methy/processed_data/tumor_clinvar.rdata")

################################# scatter plot ####################
###################### lifetime alcohol ###########
spearman_test <- cor.test(tumor_clin$PI_mc, tumor_clin$lifetimealcohol_day, method = "spearman")
spearman_cor <- spearman_test$estimate
spearman_pvalue <- spearman_test$p.value

cor(tumor_clin$PI_mc, tumor_clin$lifetimealcohol_day, method = "kendall") # 

v_mc<-ggplot(tumor_clin, aes(x=PI_mc, y=lifetimealcohol_day)) + 
  geom_point()+
  geom_smooth(method = lm)+
  labs(title = "McCartney et al. 2018",
       x = "Methylation-based alchohol consumption score",
       y = "Average lifetime daily alcohol consumption")+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman Correlation:", round(spearman_cor, 3), 
                                                   "\nP-value:", signif(spearman_pvalue, 3)), 
           hjust = 1.1, vjust = 2, size = 5, color = "red")+
  theme_minimal()

spearman_test <- cor.test(tumor_clin$PI_epitob, tumor_clin$lifetimealcohol_day, method = "spearman")
spearman_cor <- spearman_test$estimate
spearman_pvalue <- spearman_test$p.value

cor(tumor_clin$PI_epitob, tumor_clin$lifetimealcohol_day, method = "kendall") # 

v_epitob<-ggplot(tumor_clin, aes(x=PI_epitob, y=lifetimealcohol_day)) + 
  geom_point()+
  geom_smooth(method = lm)+
  labs(title = "Chamberlain et al. 2022",
       x = "Methylation-based alchohol consumption score",
       y = "Average lifetime daily alcohol consumption")+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman Correlation:", round(spearman_cor, 3), 
                                                   "\nP-value:", signif(spearman_pvalue, 3)), 
           hjust = 1.1, vjust = 2, size = 5, color = "red")+
  theme_minimal()

spearman_test <- cor.test(tumor_clin$PI_liu, tumor_clin$lifetimealcohol_day, method = "spearman")
spearman_cor <- spearman_test$estimate
spearman_pvalue <- spearman_test$p.value


v_liu<-ggplot(tumor_clin, aes(x=PI_liu, y=lifetimealcohol_day)) + 
  geom_point()+
  geom_smooth(method = lm)+
  labs(title = "Liu et al. 2018",
       x = "Methylation-based alchohol consumption score",
       y = "Average lifetime daily alcohol consumption")+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman Correlation:", round(spearman_cor, 3), 
                                                   "\nP-value:", signif(spearman_pvalue, 3)), 
           hjust = 1.1, vjust = 2, size = 5, color = "red")+
  theme_minimal()

life_time<-ggarrange(v_mc, v_epitob, v_liu, ncol = 3)

###################### last year ###########

spearman_test <- cor.test(tumor_clin$PI_mc, tumor_clin$lastyralcohol_day, method = "spearman")
spearman_cor <- spearman_test$estimate
spearman_pvalue <- spearman_test$p.value

cor(tumor_clin$PI_mc, tumor_clin$lastyralcohol_day, method = "kendall") # 

v_mc<-ggplot(tumor_clin, aes(x=PI_mc, y=lastyralcohol_day)) + 
  geom_point()+
  geom_smooth(method = lm)+
  labs(title = "McCartney et al. 2018",
       x = "Methylation-based alchohol consumption score",
       y = "Average daily alcohol consumption last year")+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman Correlation:", round(spearman_cor, 3), 
                                                   "\nP-value:", signif(spearman_pvalue, 3)), 
           hjust = 1.1, vjust = 2, size = 5, color = "red")+
  theme_minimal()

spearman_test <- cor.test(tumor_clin$PI_epitob, tumor_clin$lastyralcohol_day, method = "spearman")
spearman_cor <- spearman_test$estimate
spearman_pvalue <- spearman_test$p.value

cor(tumor_clin$PI_epitob, tumor_clin$lastyralcohol_day, method = "kendall") # 

v_epitob<-ggplot(tumor_clin, aes(x=PI_epitob, y=lastyralcohol_day)) + 
  geom_point()+
  geom_smooth(method = lm)+
  labs(title = "Chamberlain et al. 2022",
       x = "Methylation-based alchohol consumption score",
       y = "Average daily alcohol consumption last year")+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman Correlation:", round(spearman_cor, 3), 
                                                   "\nP-value:", signif(spearman_pvalue, 3)), 
           hjust = 1.1, vjust = 2, size = 5, color = "red")+
  theme_minimal()

spearman_test <- cor.test(tumor_clin$PI_liu, tumor_clin$lastyralcohol_day, method = "spearman")
spearman_cor <- spearman_test$estimate
spearman_pvalue <- spearman_test$p.value

cor(tumor_clin$PI_liu, tumor_clin$lastyralcohol_day, method = "kendall") # 

v_liu<-ggplot(tumor_clin, aes(x=PI_liu, y=lastyralcohol_day)) + 
  geom_point()+
  geom_smooth(method = lm)+
  labs(title = "Liu et al. 2018",
       x = "Methylation-based alchohol consumption score",
       y = "Average daily alcohol consumption last year")+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman Correlation:", round(spearman_cor, 3), 
                                                   "\nP-value:", signif(spearman_pvalue, 3)), 
           hjust = 1.1, vjust = 2, size = 5, color = "red")+
  theme_minimal()

recent_<-ggarrange(v_mc, v_epitob, v_liu, ncol = 3)

scatter<-ggarrange(life_time, recent_, nrow = 2)

# 8*12

################################# box plot ####################
############# lifetime consumption ############
summary(tumor_clin$liftimdrink_cat)
tumor_clin$liftimdrink_cat<-factor(as.character(tumor_clin$liftimdrink_cat),
                                   levels = c('Abstainers', 'Light drinkers', 'Moderate drinkers',
                                              'Heavy drinkers'))

tumor_clin$recentdrink_cat<-factor(as.character(tumor_clin$recentdrink_cat),
                                   levels = c('Abstainers', 'Light drinkers', 'Moderate drinkers',
                                              'Heavy drinkers'))

# standardize for comparison purpose
score<-c('PI_epitob', 'PI_mc', 'PI_liu')
tumor_clin[score]<-scale(tumor_clin[score])


summary(tumor_clin$PI_mc)
v_mc<-ggplot(tumor_clin, aes(x = liftimdrink_cat, y = PI_mc, fill = liftimdrink_cat)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) + 
  ylim(-2.5, 3)+
  labs(title = "McCartney et al. 2018",
       x = "lifetime daily alcohol consumption",
       y = "Methylation score") +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("Abstainers", "Light drinkers"), # p = 0.0000010
                                        c("Abstainers", "Moderate drinkers"),
                                        c("Abstainers", "Heavy drinkers")), # 0.0002272
                     method = "wilcox.test", # does not directly support the Tukey HSD method, add it afterwards
                     p.adjust.method = "bonferroni")


summary(tumor_clin$PI_epitob)
v_epitob<-ggplot(tumor_clin, aes(x = liftimdrink_cat, y = PI_epitob, fill = liftimdrink_cat)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) +  # narrow boxplots
  ylim(-2.5, 3)+
  labs(title = "Chamberlain et al. 2022",
       x = "lifetime daily alcohol consumption",
       y = "Methylation score") +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("Abstainers", "Light drinkers"), # p = 0.0000010
                                        c("Abstainers", "Moderate drinkers"),
                                        c("Abstainers", "Heavy drinkers")), # 0.0002272
                     method = "wilcox.test", # does not directly support the Tukey HSD method, add it afterwards
                     p.adjust.method = "bonferroni")

summary(tumor_clin$PI_liu)
v_liu<-ggplot(tumor_clin, aes(x = liftimdrink_cat, y = PI_liu, fill = liftimdrink_cat)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) + # narrow boxplots
  ylim(-2.5, 3)+ 
  labs(title = "Liu et al. 2018",
       x = "lifetime daily alcohol consumption",
       y = "Methylation score") +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("Abstainers", "Light drinkers"), # p = 0.0000010
                                        c("Abstainers", "Moderate drinkers"),
                                        c("Abstainers", "Heavy drinkers")), # 0.0002272
                     method = "wilcox.test", # does not directly support the Tukey HSD method, add it afterwards
                     p.adjust.method = "bonferroni")

life_time<-ggarrange(v_mc, v_epitob, v_liu, ncol = 3)

tumor_clin<-subset(tumor_clin, !is.na(tumor_clin$recentdrink_cat))

v_mc<-ggplot(tumor_clin, aes(x = recentdrink_cat, y = PI_mc, fill = recentdrink_cat)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) + 
  ylim(-2.5, 3)+
  labs(title = "McCartney et al. 2018",
       x = "Recent daily alcohol consumption",
       y = "Methylation score") +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("Abstainers", "Light drinkers"), # p = 0.0000010
                                        c("Abstainers", "Moderate drinkers"),
                                        c("Abstainers", "Heavy drinkers")), # 0.0002272
                     method = "wilcox.test", # does not directly support the Tukey HSD method, add it afterwards
                     p.adjust.method = "bonferroni")


summary(tumor_clin$PI_epitob)
v_epitob<-ggplot(tumor_clin, aes(x = recentdrink_cat, y = PI_epitob, fill = recentdrink_cat)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) +  # narrow boxplots
  ylim(-2.5, 3)+
  labs(title = "Chamberlain et al. 2022",
       x = "Recent daily alcohol consumption",
       y = "Methylation score") +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("Abstainers", "Light drinkers"), # p = 0.0000010
                                        c("Abstainers", "Moderate drinkers"),
                                        c("Abstainers", "Heavy drinkers")), # 0.0002272
                     method = "wilcox.test", # does not directly support the Tukey HSD method, add it afterwards
                     p.adjust.method = "bonferroni")

summary(tumor_clin$PI_liu)
v_liu<-ggplot(tumor_clin, aes(x = recentdrink_cat, y = PI_liu, fill = recentdrink_cat)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) + # narrow boxplots
  ylim(-2.5, 3)+ 
  labs(title = "Liu et al. 2018",
       x = "Recent daily alcohol consumption",
       y = "Methylation score") +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = list(c("Abstainers", "Light drinkers"), # p = 0.0000010
                                        c("Abstainers", "Moderate drinkers"),
                                        c("Abstainers", "Heavy drinkers")), # 0.0002272
                     method = "wilcox.test", # does not directly support the Tukey HSD method, add it afterwards
                     p.adjust.method = "bonferroni")

recent_1y<-ggarrange(v_mc, v_epitob, v_liu, ncol = 3)

box<-ggarrange(life_time, recent_1y, nrow = 2)

## 14*9


