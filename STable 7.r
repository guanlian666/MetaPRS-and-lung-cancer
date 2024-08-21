rm(list = ls())
gc()

library(survival)
library(plyr)

load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata") 
names(data)

## over_lung
y=Surv(time=data$difftime,event=data$outcome_LC==1)
cox_model=function(x){
    FML=as.formula(paste0('y~',x,'+data$PCA1+data$PCA2+data$PCA3+data$PCA4+data$PCA5+data$PCA6+data$PCA7+data$PCA8+data$PCA9+data$PCA10+data$chip+data$Age_when_attended_assessment_centre+data$Sex'))
    cox=coxph(FML,data=data)
    cox1=summary(cox)
    beta=round(cox1$coefficients[1,1],2)
    HR=round(cox1$coefficients[1,2],2)
    PValue=round(cox1$coefficients[1,5],3)
    CI5<-round(cox1$conf.int[1,3],2)
    CI95<-round(cox1$conf.int[1,4],2)
    uni_cox_model=data.frame('Characteristics'=x,
                            'HR'=HR,
                            'CI5'=CI5,
                            'CI95'=CI95,
                            'p'=PValue)
return(uni_cox_model)}

## Smoking
data$SI_PRS_persd=data$SI_PRS/sd(data$SI_PRS)
data$AI_PRS_persd=data$AI_PRS/sd(data$AI_PRS)
data$CPD_PRS_persd=data$CPD_PRS/sd(data$CPD_PRS)
data$SC_PRS_persd=data$SC_PRS/sd(data$SC_PRS)
data$NMR_PRS_persd=data$NMR_PRS/sd(data$NMR_PRS)
variable.names=colnames(data)[c(84:88)]
re=lapply(variable.names,cox_model)	
re2.sm=ldply(re,data.frame)   
rowname=rep('Smoking',5)
re2.sm=rbind(rowname,re2.sm)


## Lung function
data$FEV1_PRS_persd=data$FEV1_PRS/sd(data$FEV1_PRS)
data$FEV1FVC_PRS_persd=data$FEV1FVC_PRS/sd(data$FEV1FVC_PRS)
data$FVC_PRS_persd=data$FVC_PRS/sd(data$FVC_PRS)
data$PEF_PRS_persd=data$PEF_PRS/sd(data$PEF_PRS)
variable.names=colnames(data)[c(89:92)]
re=lapply(variable.names,cox_model)	
re2.LF=ldply(re,data.frame)   
rowname=rep('Lung function',5)
re2.LF=rbind(rowname,re2.LF)


##Lung diseases
data$COPD_PRS_persd=data$COPD_PRS/sd(data$COPD_PRS)
data$IPF_PRS_persd=data$IPF_PRS/sd(data$IPF_PRS)
data$ILD_PRS_persd=data$ILD_PRS/sd(data$ILD_PRS)
data$ASTHMA_PRS_persd=data$ASTHMA_PRS/sd(data$ASTHMA_PRS)
variable.names=colnames(data)[c(93:96)]
re=lapply(variable.names,cox_model)	
re2.LD=ldply(re,data.frame)   
rowname=rep('Lung function',5)
re2.LD=rbind(rowname,re2.LD)


## others
data$HEIGHT_PRS_persd=data$HEIGHT_PRS/sd(data$HEIGHT_PRS)
data$BMI_PRS_persd=data$BMI_PRS/sd(data$BMI_PRS)
data$EDU_PRS_persd=data$EDU_PRS/sd(data$EDU_PRS)
data$LCHIST_PRS_persd=data$LCHIST_PRS/sd(data$LCHIST_PRS)
variable.names=colnames(data)[c(97:100)]
re=lapply(variable.names,cox_model)	
re2.ot=ldply(re,data.frame)   
rowname=rep('others',5)
re2.ot=rbind(rowname,re2.ot)

resul=rbind(re2.sm,re2.LF,re2.LD,re2.ot)
write.csv(resul,'STable7.csv')



