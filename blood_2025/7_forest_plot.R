library(forestplot)
library(readxl)

## forest plot for self-reported BMI
selfdf <- read_excel("methy_BMI/coxforest_selfreport.xlsx")

colnames(selfdf)

ds_data<-selfdf[, c(3:5)]
ds_data$HR<-as.numeric(ds_data$HR)
ds_data$LL<-round(as.numeric(ds_data$LL), digits = 2)
ds_data$UL<-round(as.numeric(ds_data$UL), digits = 2)

colnames(ds_data)<-c("mean", "lower", "upper")
ds_data<-ds_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(selfdf))
m1fp_text<-subset(selfdf, select = c("Measurement", "HRCI"))

m1fp_table <- cbind(
  c("Measurement",  m1fp_text$Measurement),
  c("HR (95%CI)", m1fp_text$HRCI))

p1<-forestplot(m1fp_table,
               ds_data, 
               is.summary = c(TRUE, TRUE, TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,3), 
                              TRUE, TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,3),
                              TRUE, TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,3)),
               boxsize = 0.4,
               align=c("l","c","l","c","c", "c","c"),
               graph.pos = 2,
               line.margin = unit(0.6, "cm"),
               colgap = unit(2, "mm"),
               col=fpColors(box="#D73027",
                            line="#4575B4"),
               hrzl_lines = list("2" = gpar(lty = 1)),
               clip= c(0.1, 2),
               xlab = "Hazard ratio",
               zero=1,
               graphwidth=  unit(5, "cm"),
               new_page = TRUE,
               lwd.ci=1,  
               lwd.xaxis=1,  
               lwd.zero=0.5,
               vertices = TRUE,
               ci.vertices.height= 0.2,
               txt_gp=fpTxtGp(label=gpar(fontsize=12),
                              ticks=gpar(fontsize=12),
                              xlab=gpar(fontsize=12)))

# 7*9

############################################################



osdf <- read_excel("methy_BMI/coxforest_2.xlsx", 
                   sheet = "OS")

colnames(osdf)

ds_data<-osdf[, c(3:5)]
ds_data$HR<-as.numeric(ds_data$HR)
ds_data$LL<-as.numeric(ds_data$LL)
ds_data$UL<-as.numeric(ds_data$UL)

colnames(ds_data)<-c("mean", "lower", "upper")
ds_data<-ds_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(osdf))
m1fp_text<-subset(osdf, select = c("Measurement", "HRCI"))

m1fp_table <- cbind(
  c("Measurement",  m1fp_text$Measurement),
  c("HR (95%CI)", m1fp_text$HRCI))

p1<-forestplot(m1fp_table,
               ds_data, 
               is.summary = c(TRUE, TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5)),
               boxsize = 0.4,
               align=c("l","c","l","c","c", "c","c"),
               graph.pos = 2,
               line.margin = unit(0.6, "cm"),
               colgap = unit(2, "mm"),
               col=fpColors(box="#D73027",
                            line="#4575B4"),
               hrzl_lines = list("2" = gpar(lty = 1)),
               clip= c(0.1, 2),
               xlab = "Hazard ratio",
               zero=1,
               graphwidth=  unit(5, "cm"),
               new_page = TRUE,
               lwd.ci=1,  
               lwd.xaxis=1,  
               lwd.zero=0.5,
               vertices = TRUE,
               ci.vertices.height= 0.2,
               txt_gp=fpTxtGp(label=gpar(fontsize=12),
                              ticks=gpar(fontsize=12),
                              xlab=gpar(fontsize=12)))


## save 8.1 *10
cssdf <- read_excel("methy_BMI/coxforest_2.xlsx", 
                    sheet = "CSS")

colnames(cssdf)

ds_data<-cssdf[, c(3:5)]
ds_data$HR<-as.numeric(ds_data$HR)
ds_data$LL<-as.numeric(ds_data$LL)
ds_data$UL<-as.numeric(ds_data$UL)

colnames(ds_data)<-c("mean", "lower", "upper")
ds_data<-ds_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(cssdf))
m1fp_text<-subset(cssdf, select = c("Measurement", "HRCI"))

m1fp_table <- cbind(
  c("Measurement",  m1fp_text$Measurement),
  c("HR (95%CI)", m1fp_text$HRCI))

p2<-forestplot(m1fp_table,
               ds_data, 
               is.summary = c(TRUE, TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5)),
               boxsize = 0.4,
               align=c("l","c","l","c","c", "c","c"),
               graph.pos = 2,
               line.margin = unit(0.6, "cm"),
               colgap = unit(2, "mm"),
               col=fpColors(box="#D73027",
                            line="#4575B4"),
               hrzl_lines = list("2" = gpar(lty = 1)),
               clip= c(0.1, 2),
               xlab = "Hazard ratio",
               zero=1,
               graphwidth=  unit(5, "cm"),
               new_page = TRUE,
               lwd.ci=1,  
               lwd.xaxis=1,  
               lwd.zero=0.5,
               vertices = TRUE,
               ci.vertices.height= 0.2,
               txt_gp=fpTxtGp(label=gpar(fontsize=12),
                              ticks=gpar(fontsize=12),
                              xlab=gpar(fontsize=12)))

## 8.1 *10

## save 8.1 *10
cusdf <- read_excel("methy_BMI/coxforest_2.xlsx", 
                    sheet = "CUS")

colnames(cusdf)

ds_data<-cusdf[, c(3:5)]
ds_data$HR<-as.numeric(ds_data$HR)
ds_data$LL<-as.numeric(ds_data$LL)
ds_data$UL<-as.numeric(ds_data$UL)

colnames(ds_data)<-c("mean", "lower", "upper")
ds_data<-ds_data %>% add_row(mean = NA, lower = NA, upper = NA, .before = 1)

dput(names(cusdf))
m1fp_text<-subset(cusdf, select = c("Measurement", "HRCI"))

m1fp_table <- cbind(
  c("Measurement",  m1fp_text$Measurement),
  c("HR (95%CI)", m1fp_text$HRCI))

p3<-forestplot(m1fp_table,
               ds_data, 
               is.summary = c(TRUE, TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5), 
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5),
                              TRUE, rep(FALSE,5)),
               boxsize = 0.4,
               align=c("l","c","l","c","c", "c","c"),
               graph.pos = 2,
               line.margin = unit(0.6, "cm"),
               colgap = unit(2, "mm"),
               col=fpColors(box="#D73027",
                            line="#4575B4"),
               hrzl_lines = list("2" = gpar(lty = 1)),
               clip= c(0.1, 2),
               xlab = "Hazard ratio",
               zero=1,
               graphwidth=  unit(5, "cm"),
               new_page = TRUE,
               lwd.ci=1,  
               lwd.xaxis=1,  
               lwd.zero=0.5,
               vertices = TRUE,
               ci.vertices.height= 0.2,
               txt_gp=fpTxtGp(label=gpar(fontsize=12),
                              ticks=gpar(fontsize=12),
                              xlab=gpar(fontsize=12)))




