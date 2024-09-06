# The steps for calculating metaPRS and the results for Supplementary Figure 1 and Supplementary Table 8 are from this section of code
rm(list=ls())
gc()

library(data.table)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
library(dplyr)
library(glmnet)
library(pROC)

data<-as.data.frame(fread("/Public/mzm/metaPRS/data/final.data/UKB_baseline.csv"))
data_PLCO<-as.data.frame(fread("/Public/mzm/metaPRS/data/final.data/PLCO_baseline.csv"))


################################### Elastic Network Filtering in the uKB
# Model: Outcome + survival time ~21 prs+ Top 10 PCA+CHIP+AGE+SEX
# names(data[,c(58:74,76:78,81)])
# [1] "AI_PRS"      "CPD_PRS"     "NMR_PRS"     "SC_PRS"      "SI_PRS"      "FEV1_PRS"    "FVC_PRS"     "FEV1FVC_PRS" "PEF_PRS"     "ASTHMA_PRS"  "COPD_PRS"    "ILD_PRS"    
# [13] "IPF_PRS"     "BMI_PRS"     "EDU_PRS"     "LCHIST_PRS"  "HEIGHT_PRS"  "LUAD_PRS"    "LUSC_PRS"    "SCLC_PRS"    "LC_JIA_PRS"
# names(data[,c(10:11,14:24)])
# [1] "Sex"                                 "Age_when_attended_assessment_centre" "PCA1"                                "PCA2"                               
# [5] "PCA3"                                "PCA4"                                "PCA5"                                "PCA6"                               
# [9] "PCA7"                                "PCA8"                                "PCA9"                                "PCA10"                              
# [13] "chip"                               

options(digits = 10)
d1<-as.data.frame(data[,c("difftime","outcome_LC")])
A<-as.matrix(data[,c(58:74,76:78,81)])
x<-scale(A)
x<-as.data.frame(x)
x[,c(22:34)]<-as.data.frame(data[,c(10:11,14:24)])
x<-as.matrix(x)
time1<-d1$difftime
status1<-d1$outcome_LC
y<-Surv(time1,status1)


# Set alpha=0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, and select lambda.1se as the model, which can not only filter out some variables, but also have good prediction ability
#alpha=0.1
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.1)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.1,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')  # Calculate the AUC for comparison
auc(d1$outcome_LC, as.numeric(elastic_pre))
#cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.1,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))



#alpha=0.2
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.2)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.2,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
#cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.2,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))


#alpha=0.3
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.3)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.3,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
#cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.3,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))


#alpha=0.4
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.4)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.4,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
#cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.4,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))



#alpha=0.5
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.5)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
# cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.5,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
#cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.5,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))



#alpha=0.6
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.6)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.6,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
#cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.6,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))


#alpha=0.7
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.7)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
# cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.7,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
# cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.7,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))


#alpha=0.8
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.8)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
# cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.8,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
# cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.8,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))


#alpha=0.9
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.9)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
# cvfit$lambda.min
fit<-glmnet(x,y,family='cox',alpha=0.9,lambda=cvfit$lambda.min)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))
# cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha=0.9,lambda=cvfit$lambda.1se)
elastic_pre<-predict(fit,newx=x,type='response')
auc(d1$outcome_LC, as.numeric(elastic_pre))




#### Select the filtered model with alpha=0.1 (its AUC is the largest)
options(digits = 10)
d1<-as.data.frame(data[,c(43,38)])
A<-as.matrix(data[,c(58:74,76:78,81)])
x<-scale(A)
x<-as.data.frame(x)
x[,c(22:34)]<-as.data.frame(data[,c(10:11,14:24)])
x<-as.matrix(x)
time1<-d1$difftime
status1<-d1$outcome_LC
y<-Surv(time1,status1)

#alpha=0.1
set.seed(1)
cvfit=cv.glmnet(x,y,family='cox',alpha=0.1)
cvfit$lambda.min
cvfit$lambda.1se
fit<-glmnet(x,y,family='cox',alpha = 0.1)
par(mfrow=c(1,2)) 
plot(fit,xvar="lambda", label=TRUE,ylim=c(-0.075,0.25))  ###Supplementary Figure 1. A
plot(cvfit)            ###Supplementary Figure 1. B


#Extract the weights of the remaining variables after screening
coef=coef(cvfit$glmnet.fit,s=cvfit$lambda.1se,exact=F) 
coef     ### Supplementary Table 8


index=which(coef!=0) 
actcoef=coef[index]
smpvar=row.names(coef)[index] # Extract variables
smcf=cbind(var=smpvar,coef=actcoef)
smcf=as.data.frame(smcf[c(1:10),])

data<-as.data.frame(data)
B<-as.data.frame(data[smcf$var])
C<-c()
for (i in 1:ncol(B)){
  C<-rbind(C,sd(B[,i]))
  i<-i+1
} 
smcf[,3]<-C[,1]
smcf[,4]<- as.numeric(smcf[,2])/as.numeric(smcf[,3])
colnames(smcf)<-c("var","beta","sd","beta/sd")


# Calculating metaPRSs for two cohorts
###UKB
data$metaPRS<-0
train_dat<-data[smcf$var]
data$metaPRS<-apply(train_dat,1,function(x){sum(x*smcf[["beta/sd"]])})
summary(data$metaPRS)
write.csv(data,"UKB_baseline_metaPRS.csv")  #Intermediate file


###PLCO
data_PLCO$metaPRS<-0
train_dat<-data_PLCO[smcf$var]
data_PLCO$metaPRS<-apply(train_dat,1,function(x){sum(x*smcf[["beta/sd"]])})
summary(data_PLCO$metaPRS)
write.csv(data_PLCO,"PLCO_baseline_metaPRS.csv")  #Intermediate file


