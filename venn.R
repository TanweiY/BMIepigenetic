library(readxl)
scores <- read_excel("methy_BMI/BMI_cpgs.xlsx")

summary(as.factor(scores$Study))


library(VennDiagram)
scores<-subset(scores, CpGs != '(Intercept)')

# Assuming your dataframe is named df
# Convert the dataframe into a list of CpG vectors by Methylation_Score
cpg_list <- split(scores$CpGs, scores$Study)

# find the one cpgs in all scores
Reduce(intersect, cpg_list) # the following CpGs were included in all scores:
# "cg00574958 (infla)" "cg06690548 (SLC7A11)" "cg07573872 (inflammation)" "cg08548559 (inflamm)" "cg11202345" "cg14476101" "cg17501210 (inflammaotry)" "cg17901584" "cg26950531"

## find the cpgs that were in  3 groups
all_cpgs <- unlist(cpg_list)
cpg_counts <- table(all_cpgs)
names(cpg_counts[cpg_counts == 5])

# Hexadecimal color specification 
library(RColorBrewer)
brewer.pal(n = 6, name = "GnBu")

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library("ggplot2")
ggVennDiagram(cpg_list, label_alpha = 0, set_size = 4, label = c("count"),
              order.intersect.by= 'size', order.set.by = 'size',
              edge_size = 1)+
              scale_fill_gradient(low = "#CCEBC5", high = "#7BCCC4") 

# 7*5 






