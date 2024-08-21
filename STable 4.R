rm(list = ls())
gc()

library(survival)
library(plyr)
library(dplyr)


## UKB 
load('/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata')

y=Surv(time=data$difftime,event=data$outcome_LC==1)

uni_cox_model=function(x){
   model=as.formula(paste0('y~',x))
   name=variable.names
   cox=coxph(model,data=data)
   re=summary(cox)
   re2=re$concordance
   cindex=re2[1]
   cindex=round(cindex,3)
   HR = round(re$coef[1,2],2)
   HR.confint.lower <- round(re$conf.int[1,"lower .95"],2) 
   HR.confint.upper <- round(re$conf.int[1,"upper .95"],2)
   hr.ci=paste0(HR,"(", HR.confint.lower ,"-",HR.confint.upper,")")
   p.value<-signif(re$coefficients[1,"Pr(>|z|)"],4)
   p.value1=round(p.value,3)
   result=data.frame(cbind(hr.ci, p.value1,cindex))      
   return(result)}
   
names(data)
variable.names=colnames(data)[c(68:76)]
cindex=lapply(variable.names,uni_cox_model)	
cindex2=ldply(cindex,data.frame)   
rownames(cindex2)= variable.names
cindex2
 

 

## PLCO
load("/Public/mzm/metaPRS/data/final.data/PLCO_baseline_20231014.Rdata")  

y=Surv(time=data_PLCO$lung_exitdays,event=data_PLCO$lung_cancer==1)
uni_cox_model=function(x){
   model=as.formula(paste0('y~',x))
   name=variable.names
   cox=coxph(model,data=data_PLCO)
   re=summary(cox)
   re2=re$concordance
   cindex=re2[1]
   cindex=round(cindex,3)
   HR = round(re$coef[1,2],2)
   HR.confint.lower <- round(re$conf.int[1,"lower .95"],2) 
   HR.confint.upper <- round(re$conf.int[1,"upper .95"],2)
   hr.ci=paste0(HR,"(", HR.confint.lower ,"-",HR.confint.upper,")")
   p.value<-signif(re$coefficients[1,"Pr(>|z|)"],4)
   p.value1=round(p.value,3)
   result=data.frame(cbind(hr.ci, p.value1,cindex))   
   return(result)}
   
names(data_PLCO)
variable.names=colnames(data_PLCO)[c(25:33)]
cindex=lapply(variable.names,uni_cox_model)	
cindex2_1=ldply(cindex,data.frame)   
rownames(cindex2_1)= variable.names
cindex2_1

 
re=cbind(cindex2,cindex2_1)
write.csv(re,file='/Public/mzm/metaPRS/data/result/STable4.csv')

