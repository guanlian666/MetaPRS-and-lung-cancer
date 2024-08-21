rm(list = ls())
gc()
library(magrittr)

#### plco
setwd('/Public/mzm/metaPRS/data/final.data')
load('PLCO_baseline_20231014.Rdata')
load('plco.risk.lungcancer_in.plco.rdata')   ###This data comes from the intermediate file of the code in Figure 5


data_PLCO2=merge(data_PLCO,risk.lungcancer,by='#FID',all.x=T) 
data_PLCO <-data_PLCO2 

data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
model1<-glm(lung_cancer~metaPRS2+logP,family = binomial(link="logit"),data=data_PLCO)   # generate the old model logistic regression
summary(model1)

data_PLCO$prop <- predict(model1, newdata =data_PLCO, type="response")  
summary(data_PLCO$prop) 
summary(data_PLCO$P) 

  
data_PLCO$P.group=ifelse(data_PLCO$P>=0.0151,1,0)  
data_PLCO$prop.group=ifelse(data_PLCO$prop>=0.0151,1,0)  

a=with(data_PLCO,table(P.group,lung_cancer))
a2=data.frame(a) %>% .[c(3:4),]
a2$rnum=rowSums(a)
a2$cnum=colSums(a)

a2   #Export it to excel to sort out the results



data_PLCO$group=ifelse(data_PLCO$P.group==0 & data_PLCO$prop.group==0,0,
                ifelse(data_PLCO$P.group==1 & data_PLCO$prop.group==0,1,
				ifelse(data_PLCO$P.group==0 & data_PLCO$prop.group==1,2,
				ifelse(data_PLCO$P.group==1 & data_PLCO$prop.group==1,3,NA))))
library(dplyr)
b=with(data_PLCO,table(group,lung_cancer))
b2=data.frame(b) %>% .[c(5:8),]
b2$cnum=rowSums(b)
b2   #Export it to excel to sort out the results




#### ukb
load('UKB_baseline_20231018.RData')
load('plco.risk.lungcancer_in.ukb.rdata')   ###This data comes from the intermediate file of the code in Figure 5
data2=merge(data,risk.lungcancer,by='participantID')

data2$lung_cancer=data2$outcome_LC
data2$metaPRS2=data2$metaPRS
data2$prop <- predict(model1, newdata =data2, type="response")  
summary(data2$prop) 

summary(data2$P) 
 
  
data2$P.group=ifelse(data2$P>=0.0151,1,0)  
data2$prop.group=ifelse(data2$prop>=0.0151,1,0)  
with(data2,table(P.group,prop.group))

data2$group=ifelse(data2$P.group==0 & data2$prop.group==0,0,
              ifelse(data2$P.group==1 & data2$prop.group==0,1,
				ifelse(data2$P.group==0 & data2$prop.group==1,2,
				ifelse(data2$P.group==1 & data2$prop.group==1,3,NA))))

data2$lung_exitage=data2$Age_when_attended_assessment_centre+data2$difftime/365.25
a=with(data2,table(P.group,lung_cancer))
a2=data.frame(a) %>% .[c(3:4),]
a2$rnum=rowSums(a);a2   #Export it to excel to sort out the results


b=with(data2,table(group,lung_cancer))
b2=data.frame(b) %>% .[c(5:8),]
b2$cnum=rowSums(b);
b2    #Export it to excel to sort out the results





