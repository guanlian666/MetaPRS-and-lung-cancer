rm(list = ls())
gc()
library(survival)
library(plyr)
setwd('/Public/mzm/metaPRS/data')
load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata")   #data


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

##x.sd=x/sd(x)
names(data)
data$LUAD_PRS_persd=data$LUAD_PRS/sd(data$LUAD_PRS)
data$LUSC_PRS_persd=data$LUSC_PRS/sd(data$LUSC_PRS)
data$SCLC_PRS_persd=data$SCLC_PRS/sd(data$SCLC_PRS)
variable.names=colnames(data)[c(84:86)]
re=lapply(variable.names,cox_model)	
re2=ldply(re,data.frame)   
rowname=rep('Lung cancer',5)
re2=rbind(rowname,re2)

write.csv(re2,'STable5.csv')

