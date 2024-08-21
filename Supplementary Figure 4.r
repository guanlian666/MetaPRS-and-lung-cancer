rm(list = ls())
gc()

library(survival)
library(ggplot2)
library(ggsci)

setwd('/Public/mzm/metaPRS/data/')
load('PLCO_baseline_20231014.Rdata')


data_PLCO$is.smoke=ifelse(data_PLCO$cig_stat==1|data_PLCO$cig_stat==2,1,data_PLCO$cig_stat)   #0 no smoking, 1 Smoking
data_PLCO$Smoking_cessation=ifelse(data_PLCO$cig_stat==2,1,
                                   ifelse(data_PLCO$cig_stat==1,0,NA))
table(data_PLCO$Smoking_cessation,useNA = 'always')



## Deletion of missing categorical variables
da0<- subset(data_PLCO,(educat1!=9) & (COPD!=9) & (lung_fh!=9) & (race7!=99) & (ph_any_cancer!=9))  
#ph_any_cancer: Half unknown because only the intervention group filled out the dqx questionnaire
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

da0=rbind(now.smok,pre.smok,nev.smok)
risk.lungcancer=subset(da0,select=c('#FID','P','logP'))

data_PLCO2=merge(data_PLCO,risk.lungcancer,by='#FID',all.x=T) 
data_PLCO <-data_PLCO2


# ######### ######### ######### Calibration diagram
# # Generate old model logistic regression
data_PLCO3=subset(data_PLCO,is.na(data_PLCO$P)==F)
data_PLCO3$time=round((data_PLCO3$lung_exitdays/365.25),1)
con_new=data_PLCO3
con_new$lung_cancer=ifelse(con_new$time>6,0,con_new$lung_cancer)
 
Q10=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[1]]
Q20=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[2]]
Q30=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[3]]
Q40=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[4]]
Q50=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[5]]
Q60=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[6]]
Q70=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[7]]
Q80=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[8]]
Q90=quantile(con_new$P,seq(0.1,1,0.1),na.rm=T)[[9]]

con_new$AR_Q10=0
con_new[con_new$P<=Q10,]$AR_Q10=1
con_new[con_new$P>Q10&con_new$P<=Q20,]$AR_Q10=2
con_new[con_new$P>Q20&con_new$P<=Q30,]$AR_Q10=3
con_new[con_new$P>Q30&con_new$P<=Q40,]$AR_Q10=4
con_new[con_new$P>Q40&con_new$P<=Q50,]$AR_Q10=5
con_new[con_new$P>Q50&con_new$P<=Q60,]$AR_Q10=6
con_new[con_new$P>Q60&con_new$P<=Q70,]$AR_Q10=7
con_new[con_new$P>Q70&con_new$P<=Q80,]$AR_Q10=8
con_new[con_new$P>Q80&con_new$P<=Q90,]$AR_Q10=9
con_new[con_new$P>Q90,]$AR_Q10=10

results_tran <- as.data.frame(matrix(NA,ncol=4,nrow=10))
colnames(results_tran) <- c("Predicted","Observed","X95IL","X95UL")

for (i in 1:10)
{
  results_tran[i,]$Predicted=mean(con_new[con_new$AR_Q10==i,]$P)
  con_new$case=con_new$lung_cancer
  Pi<-survfit(Surv(time,case)~1, data=con_new[con_new$AR_Q10==i,],type=c("kaplan-meier"))
  results_tran[i,]$Observed=1-summary(Pi,time=10)$surv
  results_tran[i,]$X95IL=1-summary(Pi,time=10)$lower
  results_tran[i,]$X95UL=1-summary(Pi,time=10)$upper
}
 
mydata=results_tran
lm=lm(mydata$Observed~mydata$Predicted); lm 
#       intercept)  mydata$Predicted  
#      -4.423e-05         3.540e-01  

library(ResourceSelection)
a=hoslem.test(mydata$Predicted, mydata$Observed) 
 
	  
summary(lm)$r.squared
mydata$Predicted=mydata$Predicted*100
mydata$X95.IL=mydata$X95IL*100
mydata$X95.UL=mydata$X95UL*100
mydata$Observed=mydata$Observed*100


p = ggplot(mydata,aes(Predicted,ymin=X95.IL,ymax=X95.UL,color=as.factor(Predicted)))+
  xlab("Predicted probability, %")+ylab("Observed probability, %")+
  geom_point(aes(Predicted,Observed),size=3)+geom_errorbar(width=0)+scale_color_jco()+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=18,colour="black"),
        axis.title.x = element_text(size=18,colour="black"),
        axis.title.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=16,colour="black"),
        axis.text.y = element_text(size=16,colour="black"),legend.position="none",legend.title=element_blank())+
  geom_abline(intercept=-4.423e-05,slope = 3.540e-01,color="#0073C2FF",linewidth=1)+ggtitle('PLCO2014')+theme(plot.title=element_text(hjust=0.5))

ggsave(p,file='/Public/mzm/metaPRS/data/result/sFig4.pdf') 
