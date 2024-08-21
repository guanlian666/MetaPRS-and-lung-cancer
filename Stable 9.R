rm(list = ls())
gc()

load('PLCO_baseline_20231014.Rdata')   #data_PLCO
names(data_PLCO)

table(data_PLCO$cig_stop1,useNA = "always")
data_PLCO$cig_stop2=ifelse(is.na(data_PLCO$cig_stop1)==T&data_PLCO$cig_stat==0,10,data_PLCO$cig_stop1)
table(data_PLCO$cig_stop2,useNA = "always")


table(data_PLCO$arm)  #1 intervention group, 2 control group
data_PLCO$arm1=ifelse(data_PLCO$arm==2,0,data_PLCO$arm)  #0 control group，1 intervention group


###metaPRS group: Divided into 3 groups according to the top 20% and bottom 80% of metaPRS of UKB
load('UKB_baseline_20231018.RData')  #data
quq=quantile(data$metaPRS,seq(0.2,1,0.2))
quq
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
data_PLCO$PRS3=as.factor(ifelse(data_PLCO$metaPRS2<=quq[1],0,ifelse(data_PLCO$metaPRS2>quq[1] & data_PLCO$metaPRS2<quq[4],1,2)))
table(data_PLCO$PRS3)


### Calculate survival time
#Incidence of lung cancer
summary(data_PLCO$lung_exitdays)   
data_PLCO$lung_exityears=as.numeric(data_PLCO$lung_exitdays)/365.25 
data_PLCO$in_survtime=data_PLCO$lung_exityears
summary(data_PLCO$in_survtime)


#Lung cancer death
summary(data_PLCO$mortality_exitdays)
data_PLCO$mortality_exityears=as.numeric(data_PLCO$mortality_exitdays/365.25)
summary(data_PLCO$mortality_exityears)
data_PLCO$mor_survtime=data_PLCO$mortality_exityears


data_con=subset(data_PLCO,data_PLCO$arm==2)
data_inv=subset(data_PLCO,data_PLCO$arm==1)


### （1）Counting lung cancer cases and deaths
#Incidence of lung cancer
## Control group
a=table(data_con$PRS3,data_con$lung_cancer)
b=table(data_con$PRS3)
n1=data.frame(cbind(a,b))
n1$n_case=paste(n1$b,'(',n1$X1,')')

## Intervention group
c=table(data_inv$PRS3,data_inv$lung_cancer)
b=table(data_inv$PRS3)
n2=data.frame(cbind(c,b))
n2$n_case=paste(n2$b,'(',n2$X1,')')

n_1=rbind(n1,n2)


#Lung cancer death
## Control group
a=table(data_con$PRS3,data_con$f_dthl)
b=table(data_con$PRS3)
n1=data.frame(cbind(a,b))
n1$n_case=paste(n1$b,'(',n1$X1,')')

## Intervention group
c=table(data_inv$PRS3,data_inv$f_dthl)
b=table(data_inv$PRS3)
n2=data.frame(cbind(c,b))
n2$n_case=paste(n2$b,'(',n2$X1,')')

n_2=rbind(n1,n2)


library(berryFunctions)
n=rbind(n_1,n_2)
n <- insertRows(n,c(1,5,9,13),' ')
rownames(n)[1]='Control'
rownames(n)[5]='Intervention'
rownames(n)[9]='Control1'
rownames(n)[13]='Intervention1'




### （2）Calculating person year
#Incidence of lung cancer
## Control group
library(survival)
ins_con_pyears0=round(pyears(Surv(in_survtime,lung_cancer)~1,data=data_con[data_con$PRS3==0,],scale=1)$pyears,0)   ##Person-years in the low-risk group
ins_con_pyears1=round(pyears(Surv(in_survtime,lung_cancer)~1,data=data_con[data_con$PRS3==1,],scale=1)$pyears,0)   # Medium risk
ins_con_pyears2=round(pyears(Surv(in_survtime,lung_cancer)~1,data=data_con[data_con$PRS3==2,],scale=1)$pyears,0)   # High risk

## Intervention group
ins_inv_pyears0=round(pyears(Surv(in_survtime,lung_cancer)~1,data=data_inv[data_inv$PRS3==0,],scale=1)$pyears,0)  
ins_inv_pyears1=round(pyears(Surv(in_survtime,lung_cancer)~1,data=data_inv[data_inv$PRS3==1,],scale=1)$pyears,0)  
ins_inv_pyears2=round(pyears(Surv(in_survtime,lung_cancer)~1,data=data_inv[data_inv$PRS3==2,],scale=1)$pyears,0)  


#Lung cancer death
## Control group
mor_con_pyears0=round(pyears(Surv(mor_survtime,lung_cancer)~1,data=data_con[data_con$PRS3==0,],scale=1)$pyears,0)   
mor_con_pyears1=round(pyears(Surv(mor_survtime,lung_cancer)~1,data=data_con[data_con$PRS3==1,],scale=1)$pyears,0) 
mor_con_pyears2=round(pyears(Surv(mor_survtime,lung_cancer)~1,data=data_con[data_con$PRS3==2,],scale=1)$pyears,0)  

## Intervention group
mor_inv_pyears0=round(pyears(Surv(mor_survtime,lung_cancer)~1,data=data_inv[data_inv$PRS3==0,],scale=1)$pyears,0)   
mor_inv_pyears1=round(pyears(Surv(mor_survtime,lung_cancer)~1,data=data_inv[data_inv$PRS3==1,],scale=1)$pyears,0) 
mor_inv_pyears2=round(pyears(Surv(mor_survtime,lung_cancer)~1,data=data_inv[data_inv$PRS3==2,],scale=1)$pyears,0)  



### (3) Calculate the Cumulative lung cancer incidence and mortality (per 10,000 person-years)
#Incidence of lung cancer
## Control group
a=table(data_con$lung_cancer,data_con$PRS3) 
ins_con_rate0=a[2,1]/ins_con_pyears0*10000   
ins_con_rate1=a[2,2]/ins_con_pyears1*10000
ins_con_rate2=a[2,3]/ins_con_pyears2*10000
ins_con_rate=rbind(ins_con_rate0,ins_con_rate1,ins_con_rate2)

## Intervention group
b=table(data_inv$lung_cancer,data_inv$PRS3) 
ins_inv_rate0=b[2,1]/ins_inv_pyears0*10000   
ins_inv_rate1=b[2,2]/ins_inv_pyears1*10000
ins_inv_rate2=b[2,3]/ins_inv_pyears2*10000
ins_inv_rate=rbind(ins_inv_rate0,ins_inv_rate1,ins_inv_rate2)


#Lung cancer death
## Control group：
c=table(data_con$f_dthl,data_con$PRS3)  
mor_con_rate0=c[2,1]/mor_con_pyears0*10000   
mor_con_rate1=c[2,2]/mor_con_pyears1*10000
mor_con_rate2=c[2,3]/mor_con_pyears2*10000
mor_con_rate=rbind(mor_con_rate0,mor_con_rate1,mor_con_rate2)

## Intervention group
d=table(data_inv$f_dthl,data_inv$PRS3) 
mor_inv_rate0=d[2,1]/mor_inv_pyears0*10000   
mor_inv_rate1=d[2,2]/mor_inv_pyears1*10000
mor_inv_rate2=d[2,3]/mor_inv_pyears2*10000
mor_inv_rate=rbind(mor_inv_rate0,mor_inv_rate1,mor_inv_rate2)


rate=data.frame(rbind(ins_con_rate,ins_inv_rate,mor_con_rate,mor_inv_rate))
colnames(rate)=c('in_mor_rate')
rate$in_mor_rate=round(rate$in_mor_rate,2)

rate <- insertRows(rate,c(1,5,9,13),' ')
n_rate=cbind(n,rate)
n_rate=n_rate[,-c(1:3)]




####### (4) Using the low-risk group as the reference group, calculate the OR values for the other two groups
#Incidence of lung cancer
#Univariate model
res1=data.frame()
for (i in 0:1) {
  x=glm(factor(lung_cancer)~factor(PRS3),data=data_PLCO[data_PLCO$arm1==i,],family = binomial)
  OR=round(exp(x$coefficients),2)
  CI=round(exp(confint(x,level = 0.95)),2)
  OR_CI=data.frame(cbind(OR,CI))
  colnames(OR_CI)=c('OR','lower_CI','uper_CI')
  OR_CI$ORCI=paste(OR_CI[,1],'(',OR_CI$lower_CI,',',OR_CI$uper_CI,')')
  p=print(summary(x)$coefficients[,'Pr(>|z|)'],digits = 10)
  OR_CI_P=cbind(OR_CI,p)
  res1=rbind(res1,OR_CI_P)
}
res1


## Multi-variable model
res2=data.frame()
for (i in 0:1) {
  x=glm(factor(lung_cancer)~factor(PRS3)+age+factor(race7)+factor(educat1)+BMI+
          factor(COPD)+factor(ph_any_cancer)+factor(fh_cancer)+factor(cig_stat)+
          cig_stop2+cigpd_f1+cig_years1,data=data_PLCO[data_PLCO$arm1==i,],family = binomial)
  OR=round(exp(x$coefficients),2)[1:3]
  CI=round(exp(confint(x,level = 0.95)),2)[1:3,]
  OR_CI=data.frame(cbind(OR,CI))
  colnames(OR_CI)=c('OR','lower_CI','uper_CI')
  OR_CI$ORCI=paste(OR_CI[,1],'(',OR_CI$lower_CI,',',OR_CI$uper_CI,')')
  p=print(summary(x)$coefficients[,'Pr(>|z|)'],digits = 10)[1:3]
  OR_CI_P=cbind(OR_CI,p)
  res2=rbind(res2,OR_CI_P)
}
res2




#Lung cancer death
#Univariate model
res3=data.frame()
for (i in 0:1) {
  x=glm(factor(f_dthl)~factor(PRS3),data=data_PLCO[data_PLCO$arm1==i,],family = binomial)
  OR=round(exp(x$coefficients),2)
  CI=round(exp(confint(x,level = 0.95)),2)
  OR_CI=data.frame(cbind(OR,CI))
  colnames(OR_CI)=c('OR','lower_CI','uper_CI')
  OR_CI$ORCI=paste(OR_CI[,1],'(',OR_CI$lower_CI,',',OR_CI$uper_CI,')')
  p=print(summary(x)$coefficients[,'Pr(>|z|)'],digits = 10)
  OR_CI_P=cbind(OR_CI,p)
  res3=rbind(res3,OR_CI_P)
}
res3


## Multi-variable model
res4=data.frame()
for (i in 0:1) {
  x=glm(factor(f_dthl)~factor(PRS3)+age+factor(race7)+factor(educat1)+BMI+
          factor(COPD)+factor(ph_any_cancer)+factor(fh_cancer)+factor(cig_stat)+
          cig_stop2+cigpd_f1+cig_years1,data=data_PLCO[data_PLCO$arm1==i,],family = binomial)
  OR=round(exp(x$coefficients),2)[1:3]
  CI=round(exp(confint(x,level = 0.95)),2)[1:3,]
  OR_CI=data.frame(cbind(OR,CI))
  colnames(OR_CI)=c('OR','lower_CI','uper_CI')
  OR_CI$ORCI=paste(OR_CI[,1],'(',OR_CI$lower_CI,',',OR_CI$uper_CI,')')
  p=print(summary(x)$coefficients[,'Pr(>|z|)'],digits = 10)[1:3]
  OR_CI_P=cbind(OR_CI,p)
  res4=rbind(res4,OR_CI_P)
}
res4

res=rbind(cbind(res1,res2),cbind(res3,res4))
res=res[,-c(1:3,6:8)]
res[c(1,4,7,10),]=' Ref.'
res <- insertRows(res,c(1,5,9,13),' ')

res_n_rate=cbind(n_rate,res)


write.csv(res_n_rate,file = 'data_PLCO_n_rate_ins_mor_PRS3_LC_20231106.csv')


