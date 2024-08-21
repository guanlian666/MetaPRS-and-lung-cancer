
#---------------PLCO
rm(list = ls())
gc()

setwd('/Public/mzm/metaPRS/data/final.data')
load('PLCO_baseline_20231014.Rdata')
load('plco.risk.lungcancer_in.plco.rdata')

#### Data merging
data_PLCO2=merge(data_PLCO,risk.lungcancer,by='#FID',all.x=T)  
data_PLCO <-data_PLCO2 
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)


######(1) Calculate the AUC and its 95%CI
data_PLCO$case=data_PLCO$lung_cancer
model1<-glm(case~metaPRS2+logP,family = binomial(link="logit"),data=data_PLCO)  # generate the old logistic regression model
summary(model1)

library(pROC)
data_PLCO$time=round(as.numeric(data_PLCO$lung_exitdays)/365.25,2)
summary(data_PLCO$time)
data_PLCO[data_PLCO$time>6,]$case=0

data_PLCO$prop=predict(model1,data_PLCO, type = "response")
ro_0<-roc(data_PLCO$case, as.numeric(data_PLCO$logP)) 
auc=ro_0$auc   ###C-index of the PLCOm2014 in PLCO cohort
ci=ci(ro_0)    ### 95%CI for the PLCOm2014 C-index in PLCO cohort

ro_1<-roc(data_PLCO$case, as.numeric(data_PLCO$prop)) 
auc=ro_1$auc   ###C-index of the PLCOm2014+metaPRS in PLCO cohort
ci(ro_1); ci   ### 95%CI for the PLCOm2014+metaPRS C-index in PLCO cohort


#(2) Calculate categorical NRI with the nricens package
library(nricens)
table(is.na(data_PLCO$logP))
data_PLCO1<-data_PLCO[!is.na(data_PLCO$logP),]
event = ifelse(data_PLCO1$case == 1, 1, 0)

# Two matrices consisting only of predictor variables
z.std = as.matrix(subset(data_PLCO1, select = c(logP)))
z.new = as.matrix(subset(data_PLCO1, select = c(logP,metaPRS2)))

# Build two models
mo=glm(event ~., family = binomial(logit),data.frame(event,z.std),x=TRUE)
mo1=glm(event ~., family = binomial(logit),data.frame(event,z.new), x=TRUE)
summary(mo1)
p.std = mo$fitted.values
p.new = mo1$fitted.values

res=nribin(event = event, z.std = z.std, z.new = z.new,cut = 0.0151,niter =500, updown = 'category')
res    ####the categorical NRI and its 95%CI in PLCO cohort


####(3) Use survIDINRI package to calculate a continuous NRI
library(survC1)
library(survIDINRI)
table(is.na(data_PLCO$logP))
data_PLCO1<-data_PLCO[!is.na(data_PLCO$logP),]

## Big data requires a lot of running memory, and converting to matrices saves memory
dat4=as.matrix(subset(data_PLCO1,select=c("time","case","logP","metaPRS2")))

### t0 is the time to be analyzed (cannot exceed the follow-up time in the data), npert= number of iterations, covs0= old model, covs1= new model
IDI_INF<-IDI.INF(indata=dat4[,c("time","case")],
                 covs0=dat4[,c("logP")],
                 covs1=dat4[,c("logP","metaPRS2")],
                 t0=6.00,
                 npert=500,alpha=0.05)
IDI.INF.OUT(IDI_INF)    #the continuous NRI and its 95%CI in PLCO cohort





#####------------------ UKB
rm(list=ls())
gc()
#setwd('/Public/mzm/metaPRS/data/')
load('UKB_baseline_20231018.Rdata')
load('plco.risk.lungcancer_in.ukb.rdata')
data2=merge(data,risk.lungcancer,by='participantID')
data2$metaPRS2=data2$metaPRS


######(1) Calculate the AUC and its 95%CI
### Model in PLCO and predict UKB with the PLCO model
load('PLCO_baseline_20231014.Rdata')
load('plco.risk.lungcancer_in.plco.rdata')
data_PLCO2=merge(data_PLCO,risk.lungcancer,by='#FID',all.x=T) #108665
data_PLCO <-data_PLCO2 
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
data_PLCO$case=data_PLCO$lung_cancer
model1<-glm(case~metaPRS2+logP,family = binomial(link="logit"),data=data_PLCO)
summary(model1)

data2$case=data2$outcome_LC
data2$time=round(as.numeric(data2$LC_endpoint-data2$study_date)/365.25,2)
data2$case=ifelse(data2$time>6&data2$case==1,0,data2$case)
data2$prop=predict(model1,data2, type = "response")
summary(data2$prop)

library(pROC)
ro_0<-roc(data2$case, as.numeric(data2$logP)) 
auc=ro_0$auc    ### C-index of the PLCOm2014 in UKB cohort
ci(ro_0)    ### 95%CI for the PLCOm2014 C-index in UKB cohort

ro_1<-roc(data2$case, as.numeric(data2$prop)) 
auc=ro_1$auc   ###C-index of the PLCOm2014+metaPRS in ukb cohort
ci(ro_1)       ### 95%CI for the PLCOm2014+metaPRS C-index in ukb cohort
 


#(2) Calculate categorical NRI with the nricens package
rm(data_PLCO)
rm(data_PLCO2)
gc()

library(nricens)
table(is.na(data2$logP))  
data21<-data2[!is.na(data2$logP),]
event = ifelse(data21$case == 1, 1, 0)

# Two matrices consisting only of predictor variables
z.std = as.matrix(subset(data21, select = c(logP)))
z.new = as.matrix(subset(data21, select = c(logP,metaPRS2)))

# Build two models
mo=glm(event ~., family = binomial(logit),data.frame(event,z.std),x=TRUE)
mo1=glm(event ~., family = binomial(logit),data.frame(event,z.new), x=TRUE)
summary(mo1)
p.std = mo$fitted.values
p.new = mo1$fitted.values

res=nribin(event = event, z.std = z.std, z.new = z.new,cut = 0.0151,niter =500, updown = 'category')  
res    ####the categorical NRI and its 95%CI in ukb cohort



####(3) Use survIDINRI package to calculate a continuous NRI
library(survC1)
library(survIDINRI)

data21=data2[is.na(data2$logP)==F,]  
dat4=as.matrix(subset(data21,select=c("time","case","logP","metaPRS2")))

IDI_INF<-IDI.INF(indata=dat4[,c("time","case")],
                 covs0=dat4[,c("logP")],
                 covs1=dat4[,c("logP","metaPRS2")],
                 t0=6.00,
                 npert=500,alpha=0.05)
IDI.INF.OUT(IDI_INF)   #the continuous NRI and its 95%CI in PLCO cohort






