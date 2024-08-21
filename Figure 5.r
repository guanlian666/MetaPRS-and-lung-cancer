rm(list = ls())
gc()

setwd('/Public/mzm/metaPRS/data/final.data')
load('PLCO_baseline_20231014.Rdata')
names(data_PLCO)


table(data_PLCO$cig_stat,useNA = 'always')
data_PLCO$is.smoke=ifelse(data_PLCO$cig_stat==1|data_PLCO$cig_stat==2,1,data_PLCO$cig_stat)   #0 no smoking, 1 Smoking

data_PLCO$Smoking_cessation=ifelse(data_PLCO$cig_stat==2,1,
                                   ifelse(data_PLCO$cig_stat==1,0,NA))
table(data_PLCO$Smoking_cessation,useNA = 'always')


## Deletion of missing categorical variables
da0<- subset(data_PLCO,(educat1!=9) & (COPD!=9) & (lung_fh!=9) & (race7!=99) & (ph_any_cancer!=9))  
da0$Black=ifelse(da0$race7==2,1,0)
da0$Hispanic=ifelse(da0$race7==3,1,0)
da0$Asian=ifelse(da0$race7==4,1,0)
da0$Pacific=ifelse(da0$race7==5,1,0)
da0$Native=ifelse(da0$race7==6,1,0) 
 
dim(da0) 
dim(data_PLCO)

## smoke model
now.smok=da0[da0$cig_stat==1,]
attach(now.smok)
now.smok$logP=-6.84088+(age-62)*0.0796252-0.0879622*(educat1-4)-(BMI-27)*0.0289916+(COPD)*0.3454183+0.3213965*(Black)-0.8202554*(Hispanic)-0.5240639 *(Asian)-1.364461*(Pacific) + 0.9521109*(Native)+ 0.4845352*(ph_any_cancer)+0.5857166*(lung_fh)+2.890119*1-0.1868627*(((cigpd_f1+0.25)/100)^(-1)-4)+0.0305386*cig_years1
now.smok$P=exp(now.smok$logP)/(1+exp(now.smok$logP))
detach(now.smok)
summary(now.smok$P)


## pre-smoke model
pre.smok=da0[da0$cig_stat==2,]
attach(pre.smok)
pre.smok$logP=-6.84088+(age-62)*0.0796252-0.0879622*(educat1-4)-(BMI-27)*0.0289916+(COPD)*0.3454183+0.3213965*(Black)-0.8202554*(Hispanic)-0.5240639 *(Asian)-1.364461*(Pacific) + 0.9521109*(Native)+ 0.4845352*(ph_any_cancer)+0.5857166*(lung_fh)+2.310767*1-0.1868627*(((cigpd_f1+0.25)/100)^(-1)-4)+0.0305386*cig_years1-0.0321588*(cig_stop1-10)
pre.smok$P=exp(pre.smok$logP)/(1+exp(pre.smok$logP))
detach(pre.smok)
summary(pre.smok$P)
 

## never-smoke model
nev.smok=da0[da0$cig_stat==0,]
attach(nev.smok)
nev.smok$logP=-6.84088+(age-62)*0.0796252-0.0879622*(educat1-4)-(BMI-27)*0.0289916+(COPD)*0.3454183+0.3213965*(Black)-0.8202554*(Hispanic)-0.5240639 *(Asian)-1.364461*(Pacific) + 0.9521109*(Native)+ 0.4845352*(ph_any_cancer)+ 0.5857166*(lung_fh)
nev.smok$P=exp(nev.smok$logP)/(1+exp(nev.smok$logP))
detach(nev.smok)
summary(nev.smok$P)
 

da0=rbind(now.smok,pre.smok,nev.smok) # 102905 
risk.lungcancer=subset(da0,select=c('#FID','logP','P'))
save(risk.lungcancer,file='plco.risk.lungcancer_in.plco.rdata')     


#### data merging
library(magrittr)
data_PLCO2=merge(data_PLCO,risk.lungcancer,by='#FID',all.x=T)  
data_PLCO <-data_PLCO2 


###metaPRS group: Divided into 3 layers according to the top 20% and bottom 80% of metaPRS in UKB cohort
load('UKB_baseline_20231018.Rdata')
quq=quantile(data$metaPRS,seq(0.2,1,0.2))
quq

data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)

attach(data_PLCO)
data_PLCO$group=as.factor(ifelse(metaPRS2<=quq[1],1,ifelse(metaPRS2>quq[1] & metaPRS2<quq[4],2,3)))
detach(data_PLCO)

data_PLCO$lung_exitdays=data_PLCO$lung_exitdays/365.25


###The population was grouped according to different models and genetic risk
data_PLCO$group2=ifelse(data_PLCO$group==1 & data_PLCO$P<0.0134,1,
                 ifelse(data_PLCO$group==2 & data_PLCO$P<0.0134,2,
                 ifelse(data_PLCO$group==3 & data_PLCO$P<0.0134,3,

                 ifelse(data_PLCO$group==1 & (data_PLCO$P>=0.0134 & data_PLCO$P<0.0151),4,
                 ifelse(data_PLCO$group==2 & (data_PLCO$P>=0.0134 & data_PLCO$P<0.0151),5,
                 ifelse(data_PLCO$group==3 & (data_PLCO$P>=0.0134 & data_PLCO$P<0.0151),6,

				         ifelse(data_PLCO$group==1 & (data_PLCO$P>=0.0151 & data_PLCO$P<0.15),7,
                 ifelse(data_PLCO$group==2 & (data_PLCO$P>=0.0151 & data_PLCO$P<0.15),8,
                 ifelse(data_PLCO$group==3 & (data_PLCO$P>=0.0151 & data_PLCO$P<0.15),9, 
				 
				         ifelse(data_PLCO$group==1 & (data_PLCO$P>=0.15),10,
                 ifelse(data_PLCO$group==2 & (data_PLCO$P>=0.15),11,
                 ifelse(data_PLCO$group==3 & (data_PLCO$P>=0.15),12, 
				 NA))))))))))))
				 
table(data_PLCO$group2,useNA = 'always')				 


## Predict risks
model1<-glm(lung_cancer~metaPRS+logP,family = binomial(link="logit"),data=data_PLCO)
summary(model1)
data_PLCO$time=round((data_PLCO$lung_exitdays),1) 
con_new=data_PLCO
con_new$lung_cancer=ifelse(con_new$time>6,0,con_new$lung_cancer)
con_new$prop=predict(model1,con_new, type = "response")


## Smokers
da_s=con_new[con_new$is.smoke==1,]
table(da_s$group2)  #12 groups
res1=data.frame()
for (i in 1:12) {
  mean=mean(da_s[da_s$group2==i,]$prop,na.rm=T)
  sd= sd(da_s[da_s$group2==i,]$prop, na.rm = TRUE)
  res=cbind(mean,sd)
  res1=rbind(res1,res)
}
group2=c(1:12)
res1=cbind(group2,res1)
res1


##Non smokers
da_n=con_new[con_new$is.smoke==0,] 
table(da_n$group2)   #Only three groups
res2=data.frame()
for (i in 1:3) {
  mean=mean(da_n[da_n$group2==i,]$prop,na.rm=T)
  sd= sd(da_n[da_n$group2==i,]$prop, na.rm = TRUE)
  res=cbind(mean,sd)
  res2=rbind(res2,res)
}
group2=c(1:3)
res2=cbind(group2,res2)
res2



### Plot the results for the smoking population
library(ggplot2)
library(plyr)
library(patchwork)
plco=res1[-13,]
plco$mean=round(plco$mean*100,2)
plco$sd=plco$sd*100
plco$group=c(rep(c('Low','Intermediate','High'),4))
plco$group1=rep(c('Model-low','Model-intermediate','Model-high','Model-very high'),each=3)
plco$group <-factor(plco$group,ordered=TRUE,levels=c('Low','Intermediate','High')) 


fig.1=ggplot(plco, aes(x = group1 , y = mean, fill = group)) +
  geom_bar(
    stat = "identity",
    aes(fill=group),
    position = position_dodge(),
    width = 0.7
  ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7),
                width = 0.3) +
  geom_text(aes(label=mean),size=3.5,vjust=-0.5,position =position_dodge(0.7))+ 
  geom_hline(yintercept = 1.51,colour="grey", linetype="dashed")+
  scale_fill_manual(values = c("#72c574", "#fd8c3c", "#e51a1d","#72c574", "#fd8c3c", "#e51a1d")) +
  scale_x_discrete(limits=c('Model-low','Model-intermediate','Model-high','Model-very high'),expand = c(0, 0.35))+  
  scale_y_continuous(breaks = seq(0, 25, 5), 
                     limits = c(0, 25),
                     expand = c(0, 0)) +
  labs(x = "Clinical risk category by PLCO2014 model", y = "6–year absolute risk (%)", fill = 'MetaPRS') +
  theme_classic()+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.position = c(0.2,0.8))

  
  
#####------------------ UKB
load('UKB_baseline_20231018.Rdata')
load('asthma_COPD_IPF_ILD_ICD10_method2.rdata')
data2=merge(data,data9,by='participantID',all.x=T)
data2$COPD=ifelse(data2$COPD_ICD10==1,1,0)

da0<- subset(data2,(educate_2!=9) & (COPD!=9) & (FH_Lung!=9) & (Ethnic_background2!=99))  

da0$Black=ifelse(da0$Ethnic_background2==2,1,0)
da0$Asian=ifelse(da0$Ethnic_background2==4,1,0)

## smoke model
now.smok=da0[da0$Smoking_status==2,]
attach(now.smok)
now.smok$logP=-6.84088+(Age_when_attended_assessment_centre-62)*0.0796252-0.08796400*(educate_2-4)-(BMI-27)*0.0289916+(COPD)*0.3454183+0.34003965*(Black)-0.5240639 *(Asian) +0.5857166*(FH_Lung)+2.890119*1-0.1868627*(((Number_of_cigarettes_daily_sex_ever+0.25)/100)^(-1)-4)+0.04005386*smoke_year
now.smok$P=exp(now.smok$logP)/(1+exp(now.smok$logP))
detach(now.smok)
summary(now.smok$P)


## pre-smoke model
pre.smok=da0[da0$Smoking_status==1,]
attach(pre.smok)
pre.smok$logP=-6.84088+(Age_when_attended_assessment_centre-62)*0.0796252-0.08796400*(educate_2-4)-(BMI-27)*0.0289916+(COPD)*0.3454183+0.34003965*(Black)-0.5240639 *(Asian) +0.5857166*(FH_Lung)+2.310767*1-0.1868627*(((Number_of_cigarettes_daily_sex_ever+0.25)/100)^(-1)-4)+0.04005386*smoke_year-0.03400588*(Smoking_quit_time-10)
pre.smok$P=exp(pre.smok$logP)/(1+exp(pre.smok$logP))
detach(pre.smok)
summary(pre.smok$P)


## never-smoke model
nev.smok=da0[da0$Smoking_status==0,]
attach(nev.smok)
nev.smok$logP=-6.84088+(Age_when_attended_assessment_centre-62)*0.0796252-0.08796400*(educate_2-4)-(BMI-27)*0.0289916+(COPD)*0.3454183+0.34003965*(Black)-0.5240639 *(Asian) + 0.5857166*(FH_Lung)
nev.smok$P=exp(nev.smok$logP)/(1+exp(nev.smok$logP))
detach(nev.smok)
summary(nev.smok$P)
 
da0=rbind(now.smok,pre.smok,nev.smok) # 102905 
risk.lungcancer=subset(da0,select=c('participantID','logP','P'))
save(risk.lungcancer,file='plco.risk.lungcancer_in.ukb.rdata')


data=merge(data,risk.lungcancer,by='participantID')
data$lung_exitdays=data$difftime/365.25
quq=quantile(data$metaPRS,seq(0.2,1,0.2))
quq

data$group=as.factor(ifelse(data$metaPRS<=quq[1],1,ifelse(data$metaPRS>quq[1] & data$metaPRS<quq[4],2,3)))
table(data$group,useNA = 'always')


###The population was grouped according to different models and genetic risk
data$group2=ifelse(data$group==1 & data$P<0.0134,1,
                   ifelse(data$group==2 & data$P<0.0134,2,
                          ifelse(data$group==3 & data$P<0.0134,3,
                                 
                                 ifelse(data$group==1 & (data$P>=0.0134 & data$P<0.0151),4,
                                        ifelse(data$group==2 & (data$P>=0.0134 & data$P<0.0151),5,
                                               ifelse(data$group==3 & (data$P>=0.0134 & data$P<0.0151),6,
                                                      
                                                      ifelse(data$group==1 & (data$P>=0.0151 & data$P<0.15),7,
                                                             ifelse(data$group==2 & (data$P>=0.0151 & data$P<0.15),8,
                                                                    ifelse(data$group==3 & (data$P>=0.0151 & data$P<0.15),9, 
                                                                           
                                                                           ifelse(data$group==1 & (data$P>=0.15),10,
                                                                                  ifelse(data$group==2 & (data$P>=0.15),11,
                                                                                         ifelse(data$group==3 & (data$P>=0.15),12, 
                                                                                                NA))))))))))))
table(data$group2,useNA = 'always')


## Predict risk
data_PLCO$metaPRS=data_PLCO$metaPRS-(7.453-0.157561)
model1<-glm(lung_cancer~metaPRS+logP,family = binomial(link="logit"),data=data_PLCO) # Generate old model logistic regression
summary(model1)

data$time=round(as.numeric(data$LC_endpoint-data$study_date)/365.25,1)
con_new=data
con_new$outcome_LC1=con_new$outcome_LC
con_new$outcome_LC=ifelse(con_new$time>6,0,con_new$outcome_LC)
con_new$prop=predict(model1,con_new, type = "response")


## Smoker
con_new$is.smoke=ifelse(con_new$Smoking_status!=0,1,0)
da_s=con_new[con_new$is.smoke==1,]

table(da_s$group2)  #12 groups
res1=data.frame()
for (i in 1:12) {
  mean=mean(da_s[da_s$group2==i,]$prop,na.rm=T)
  sd= sd(da_s[da_s$group2==i,]$prop, na.rm = TRUE)
  res=cbind(mean,sd)
  res1=rbind(res1,res)
}
group2=c(1:12)
res1=cbind(group2,res1)
res1

### Non-smoker
da_n=con_new[con_new$is.smoke==0,]
table(da_n$group2)  #3 groups
res2=data.frame()
for (i in 1:3) {
  mean=mean(da_n[da_n$group2==i,]$prop,na.rm=T)
  sd= sd(da_n[da_n$group2==i,]$prop, na.rm = TRUE)
  res=cbind(mean,sd)
  res2=rbind(res2,res)
}
group2=c(1:3)
res2=cbind(group2,res2)
res2



#### Plot the results for the smoking population
ukb=res1
ukb$mean=round(ukb$mean*100,2)
ukb$sd=ukb$sd*100
ukb$group=c(rep(c('Low','Intermediate','High'),4))
ukb$group1=rep(c('Model-low','Model-intermediate','Model-high','Model-very high'),each=3)
ukb$group <-factor(ukb$group,ordered=TRUE,levels=c('Low','Intermediate','High')) #修改因子水平 
ukb

b=ggplot(ukb, aes(x = group1 , y = mean, fill = group)) +
  geom_bar(
    stat = "identity",
    aes(fill=group),
    position = position_dodge(),
    #color = "white",
    width = 0.7
  ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.7),
                width = 0.3) +
  geom_text(aes(label=mean),size=3.5,vjust=-0.5,position =position_dodge(0.7))+  
  geom_hline(yintercept = 1.51,colour="grey", linetype="dashed")+
  scale_fill_manual(values = c("#72c574", "#fd8c3c", "#e51a1d","#72c574", "#fd8c3c", "#e51a1d")) +
  scale_x_discrete(limits=c('Model-low','Model-intermediate','Model-high','Model-very high'),expand = c(0, 0.35))+  
  scale_y_continuous(breaks = seq(0, 25, 5), 
                     limits = c(0, 25),
                     expand = c(0, 0)) +
  labs(x = "Clinical risk category by PLCO2014 model in UKB", y = "6–year absolute risk (%)", fill = 'MetaPRS') +
  theme_classic()+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.position = c(0.2,0.8))

b


fig=fig.1+b+plot_annotation('A')
ggsave('/Public/mzm/metaPRS/data/result/fig4.pdf')