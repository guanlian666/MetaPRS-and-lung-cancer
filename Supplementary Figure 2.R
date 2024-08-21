rm(list=ls())
gc()

setwd('/Public/mzm/metaPRS/data/final.data')
load('PLCO_baseline_20231014.Rdata')   #data_PLCO


###metaPRS grouping: Divided into 3 layers according to the top 20% and bottom 80% of metaPRS in UKB cohort
load('UKB_baseline_20231018.Rdata')
quq=quantile(data$metaPRS,seq(0.2,1,0.2))
quq
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
data_PLCO$PRS3=as.factor(ifelse(data_PLCO$metaPRS2<=quq[1],0,ifelse(data_PLCO$metaPRS2>quq[1] & data_PLCO$metaPRS2<quq[4],1,2)))




#### ---------------------(1)Incidence of lung cancer
LC_ins=subset(data_PLCO,data_PLCO$lung_cancer==1)  

library(dplyr) 
library(ggplot2)  
library(cowplot) 
library(ggpubr)  
library(gridExtra) 


# Create a scatter plot
p1 <- ggscatter(LC_ins, x = "metaPRS2", y = "lung_exitage",
               size = 1.5,
               add = "reg.line",  
               add.params = list(color = "#77C034", fill = "#C5E99B", size = 1),  # 自定义回归线的颜色
               conf.int = TRUE) +   
  stat_cor(method = "pearson", label.x = -1.35, label.y = 87, label.sep = "\n") +
  xlab("MetaPRS") +
  ylab("Age of onset (years)")




####-------------- (2) Lung cancer deaths
LC_mor=subset(data_PLCO,data_PLCO$f_dthl==1)   
p2 <- ggscatter(LC_mor, x = "metaPRS2", y = "mortality_exitage",
               size = 1.5,
               add = "reg.line",  
               add.params = list(color = "#77C034", fill = "#C5E99B", size = 1),  # 自定义回归线的颜色
               conf.int = TRUE) +   
  stat_cor(method = "pearson", label.x = -1.2, label.y = 92, label.sep = "\n") +
  xlab("MetaPRS") +
  ylab("Age at death (years)")



library(patchwork)
p1+p2+plot_annotation(tag_levels='A')
ggsave('/Public/mzm/metaPRS/data/result/SFigure2.pdf')


