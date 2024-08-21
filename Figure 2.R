rm(list=ls())
gc()
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggcorrplot")
screen -r 11175.m

library("ggplot2")
library("ggstatsplot")
library("ggcorrplot") 
library('corrplot')
 
setwd('/Public/mzm/metaPRS/data/final.data')
load('UKB_baseline_20231018.Rdata')   #data

names(data)
data_PRS<-data[,c("LC_JIA_PRS","LUAD_PRS","LUSC_PRS","SCLC_PRS",
                   "AI_PRS","SI_PRS","SC_PRS","CPD_PRS","NMR_PRS",
                   "FEV1_PRS","FVC_PRS","FEV1FVC_PRS","PEF_PRS",
                   "COPD_PRS","ASTHMA_PRS","ILD_PRS","IPF_PRS",
                   "BMI_PRS","EDU_PRS","LCHIST_PRS","HEIGHT_PRS")]
				   
colnames(data_PRS)<-c("LC-PRS","LUAD-PRS","LUSC-PRS","SCLC-PRS",
                      "AI-PRS","SI-PRS","SC-PRS","CPD-PRS","NMR-PRS",
                      "FEV1-PRS","FVC-PRS","FEV1/FVC-PRS","PEF-PRS",
                      "COPD-PRS","ASTHMA-PRS","ILD-PRS","IPF-PRS",
                      "BMI-PRS","Education-PRS","FHLC-PRS","Height-PRS") 
				   
ggcorrmat(
  data = data_PRS,
  type = "pearson",
  colors = c("blue", "white", "red"),
  matrix.type = "upper", 
  pch = "cross", #Shape of non-significant marker (when p≥0.05)
  ggcorrplot.args = list(method = "square",
                         outline.color = "black",
                         hc.order = F))



