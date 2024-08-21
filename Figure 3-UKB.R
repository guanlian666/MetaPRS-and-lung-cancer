
###########(1)figure 3-A
rm(list = ls())
gc()

library(ggplot2)
library(patchwork)

load('/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata')
data$outcome_LC=as.numeric(data$outcome_LC)
data$outcome_LC=factor(data$outcome_LC,levels=c(0,1),labels=c("No","Yes"))

pp=ggplot(data,aes(metaPRS,fill=outcome_LC))+
  geom_density(position='identity',alpha=0.8)+
  xlab('metaPRS')+scale_fill_manual(values=c("#00468B99","#AD002A99"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_line(colour='black'),
        axis.text=element_text(size=8,colour = "black"),
        axis.title=element_text(size=10,colour='black',face='bold'))+
  theme(legend.position=c(0.7,1),legend.key.size=unit(0.8,'lines'))+
  ggtitle('UKB')+theme(plot.title=element_text(hjust=0.5))
pp




###########(2)figure 3-B
rm(list = ls())
gc()
library(survival)
library(plyr)
library(dplyr)

setwd('/Public/mzm/metaPRS/data/final.data')
load('UKB_baseline_20231018.RData')
y=Surv(time=data$difftime,event=data$outcome_LC==1)
uni_cox_model=function(x){
  model=as.formula(paste0('y~',x))
  name=variable.names
  cox=coxph(model,data=data)
  re=summary(cox)
  re2=re$concordance
  cindex=re2[1]
  cindex=round(cindex,3)
  se=re2[2]
  result=data.frame(cbind(cindex,se))      
  return(result)}

variable.names=c('metaPRS',"LC_DAI_PRS", "LC_GRAFF_PRS","LC_JIA_PRS","LC_SHI_PRS","LC_ZHANG_PRS",
                 "LC_FRITSCHE_14" ,"LC_FRITSCHE_19","LC_HUNG_35","LC_HUNG_128")
cindex=lapply(variable.names,uni_cox_model)	
cindex1=ldply(cindex,data.frame)   
rownames(cindex1)= variable.names


######plot figure: Rstudio
library(ggsci)
library(gridExtra)
library(ggplot2)
library(patchwork)
c_index=cindex1
c_index$X1=c("MetaPRS","LC-DAI-PRS","LC-GRAFF-PRS","LC-JIA-PRS","LC-SHI-PRS","LC-ZHANG-PRS",
             "LC-FRITSCHE-14","LC-FRITSCHE-19","LC-HUNG-35","LC-HUNG-128")
fig.p0=ggplot(c_index, aes(x = X1, y = cindex, fill = X1, color = "black")) + 
  ylim(0,0.7)+
  geom_col(position = "dodge", color = "#1f77b4", width = 0.8) +
  geom_errorbar(aes(x = X1, ymin=cindex-se, ymax=cindex+se),width = 0.3,linewidth = 1, position = position_dodge(0.8),color = "black") + 
  labs(x = "PRS", y = "C-index", fill = "PRS") +
  scale_fill_manual(values = c("#1f77b4","#1f77b4","#1f77b4","#1f77b4","#1f77b4","#1f77b4","#1f77b4","#1f77b4","#1f77b4","#1f77b4"))+
  theme_bw()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12,face = "bold"))+
  theme(axis.text.y = element_text(size = 12,face = "bold"))+
  theme(axis.title.x = element_text(size = 14, face = "bold")) +
  theme(axis.title.y = element_text(size = 14, face = "bold"))+
  ggtitle('UKB')+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none")
fig.p0




##############(3)figure 3-C
rm(list = ls())
gc()
library(splines)
library(lattice)
library(rms)
library(ggplot2)
library(survminer)
library(survival)
library(data.table)

load('/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata')
data$metaPRS.sd=data$metaPRS/sd(data$metaPRS)
pcubic=datadist(data)
options(datadist='pcubic')
fit_3=cph(Surv(difftime,outcome_LC) ~ rcs(metaPRS.sd,3)+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+
            chip+Age_when_attended_assessment_centre+Sex, data = data,x=T,y=T)
fit_4=cph(Surv(difftime,outcome_LC) ~ rcs(metaPRS.sd,4)+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+
            chip+Age_when_attended_assessment_centre+Sex, data = data,x=T,y=T)
fit_5=cph(Surv(difftime,outcome_LC) ~ rcs(metaPRS.sd,5)+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+
            chip+Age_when_attended_assessment_centre+Sex, data = data,x=T,y=T)
AIC(fit_3) 
AIC(fit_4) 
AIC(fit_5) 


# The ref is adjusted according to y-hat =1, specifying a reference value to package the data again, OR=1
Pre_HR <-rms::Predict(fit_3,metaPRS.sd,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
ggplot(Pre_HR)
# Select the metaPRS.sd value corresponding to y-hat=1
a=subset(Pre_HR,round(Pre_HR$yhat,0)==1,select=c('yhat','metaPRS.sd')) #0.62431978

ddist <- datadist(data)
ddist$limits$metaPRS.sd[2] <- 0.62431978
options(datadist="ddist")

fit_3=cph(Surv(difftime,outcome_LC) ~ rcs(metaPRS.sd,3)+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+
            chip+Age_when_attended_assessment_centre+Sex, data = data,x=T,y=T)
anova(fit_3)   #The p-value in the figure 3-C
pred=Predict(fit_3,metaPRS.sd,fun=exp,ref.zero=T)


### figure
PP0=ggplot(pred,ylim=c(0,5))+
  geom_ribbon(aes(metaPRS.sd,ymin = lower, ymax = upper),linetype='blank',linewidth=0.7,fill="lightskyblue1")+theme_bw()+
  geom_line(aes(x=metaPRS.sd ,y=yhat),linetype=1,linewidth=1,color="#0099CC")+
  geom_hline(yintercept=1, linetype=2,color="#AD002A99",linewidth=1)+
  labs(x="metaPRS (per sd increased)", y="HR of lung cancer")+
  theme(panel.background=element_blank(),panel.grid=element_blank(),
        axis.text=element_text(size=10,colour='black'),axis.title.x=element_text(size=12,
                                                                                 vjust=-3),axis.title.y=element_text(size=12,vjust=5, angle = 90),
        plot.caption=element_blank(),plot.margin = unit(c(1, 1, 1, 1),"cm"))
PP0


# HR on the figure
output=function(x){   x <- summary(x)
HR = round(x$coef[1,2],2)
HR.confint.lower <- round(x$conf.int[1,"lower .95"],2) 
HR.confint.upper <- round(x$conf.int[1,"upper .95"],2)
p.value<-signif(x$coefficients[1,"Pr(>|z|)"],4)
p.value1=round(p.value,3)
res<-c(HR,HR.confint.lower, HR.confint.upper,p.value1)
return(res)}
fit_3=coxph(Surv(difftime,outcome_LC) ~ metaPRS.sd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+
              chip+Age_when_attended_assessment_centre+Sex, data = data)
res=output(fit_3);res





#######(4)figure 3-D
rm(list = ls())
gc()

library(ggplot2)
library(survminer)
library(survival)
library(data.table)

load('/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.RData')
data$metaPRS.sd=data$metaPRS/sd(data$metaPRS)

###Distribution: The population was divided into 3 groups according to the quartile
quq=quantile(data$metaPRS.sd,seq(0.2,1,0.2))
quq

attach(data)
data$group=as.factor(ifelse(metaPRS.sd<=quq[1],1,ifelse(metaPRS.sd>quq[1] & metaPRS.sd<quq[4],2,3)))
detach(data)

data$difftime2=data$difftime/365.25
output1=function(x){   x <- summary(x)
HR = round(x$coef[c(1:2),2],2)
HR.confint.lower <- round(x$conf.int[c(1:2),"lower .95"],2) 
HR.confint.upper <- round(x$conf.int[c(1:2),"upper .95"],2)
p.value<-signif(x$coefficients[c(1:2),"Pr(>|z|)"],3)
p.value1=round(p.value,3)
res<-cbind(HR,HR.confint.lower, HR.confint.upper,p.value1)
return(res)}

fit_3=coxph(Surv(difftime2,outcome_LC) ~ factor(group)+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+
              factor(chip)+Age_when_attended_assessment_centre+factor(Sex), data = data)
res=output1(fit_3);res   #HR and 95%CI in figure


#### Data standardization
dat4=subset(data,select=c('group', 'Sex', 'Age_when_attended_assessment_centre',
                          'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6', 'PCA7', 'PCA8', 'PCA9', 'PCA10', 'chip' ))
dat_df=with(dat4,
            data.frame(difftime2=rep(10,111),
                       outcome_LC=rep(1,111),
                       Sex=rep(round(median(Sex,na.rm=TRUE),0),111),
                       Age_when_attended_assessment_centre=rep(c(37:73),3), 
                       PCA1=rep(mean(PCA1,na.rm=TRUE),111),
                       PCA2=rep(mean(PCA2,na.rm=TRUE),111),
                       PCA3=rep(mean(PCA3,na.rm=TRUE),111),
                       PCA4=rep(mean(PCA4,na.rm=TRUE),111),
                       PCA5=rep(mean(PCA5,na.rm=TRUE),111),
                       PCA6=rep(mean(PCA6,na.rm=TRUE),111),
                       PCA7=rep(mean(PCA7,na.rm=TRUE),111),
                       PCA8=rep(mean(PCA8,na.rm=TRUE),111),
                       PCA9=rep(mean(PCA9,na.rm=TRUE),111),
                       PCA10=rep(mean(PCA10,na.rm=TRUE),111),
                       chip=rep(round(mean(chip,na.rm=TRUE),0),111),
                       group=c(rep('1',37),rep('2',37),rep('3',37))))
dat_df$group=as.factor(dat_df$group)
dat_df$Sex=as.factor(dat_df$Sex)
dat_df$chip=as.factor(dat_df$chip)


## Calculate survival time
p <- predict(fit_3, newdata = dat_df, type="expected")  

#plot Figures#
plotdata=cbind(dat_df,p)
names(plotdata)[17]=c("abosulate_rate")
plotdata$abosulate_rate=plotdata$abosulate_rate*100

p1<-ggplot(plotdata,aes(x=Age_when_attended_assessment_centre,y=abosulate_rate, colour=group,linetype=group))+
  geom_line(linewidth=0.8)+
  labs(x="Age",y="Absolute risk,%")+
  theme(panel.grid.major=element_line(colour=NA))+
  scale_linetype_manual(values=c('solid','solid','solid'))+
  scale_color_manual(values = c('#86C27D',"#4A6990FF","#DB423E"),labels=c("Low genetic risk",
                                                                          "Intermediate genetic risk: 1.68(1.50-1.87)",'High genetic risk:2.52(2.23-2.84)'))+
  theme(legend.background = element_blank(),panel.background=element_rect(fill='transparent'))+
  theme(panel.border=element_rect(color = "grey", fill = NA, linewidth = 0.8))+
  theme(legend.position = c(0.1, .95),
        legend.justification = c("left", "top"),panel.grid.major = element_line(colour = "grey", linetype = "dashed"))

p1









