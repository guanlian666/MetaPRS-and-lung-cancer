
########(1)figure 4-A
rm(list = ls())
gc()

library(ggplot2)
library(patchwork)


load("/Public/mzm/metaPRS/data/final.data/PLCO_baseline_20231014.Rdata")  
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)

data_PLCO$lung_cancer=factor(data_PLCO$lung_cancer,levels=c(0,1),labels=c("No","Yes"))

pp2=ggplot(data_PLCO,aes(metaPRS2,fill=lung_cancer))+
  geom_density(position='identity',alpha=0.8)+
  xlab('metaPRS')+scale_fill_manual(values=c("#00468B99","#AD002A99"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_line(colour='black'),
        axis.text=element_text(size=8,colour = "black"),
        axis.title=element_text(size=10,colour='black',face='bold'))+
  theme(legend.position=c(0.7,1),legend.key.size=unit(0.8,'lines'))+
  ggtitle('PLCO')+theme(plot.title=element_text(hjust=0.5))
pp2



#####(2)figure 4-B
rm(list = ls())
gc()

library(survival)
library(plyr)
library(dplyr)

load("PLCO_baseline_20231014.Rdata")  

y=Surv(time=data_PLCO$lung_exitdays,event=data_PLCO$lung_cancer==1)
uni_cox_model=function(x){
  model=as.formula(paste0('y~',x))
  name=variable.names
  cox=coxph(model,data=data_PLCO)
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
cindex2=ldply(cindex,data.frame)   
rownames(cindex2)= variable.names

 
c_index_plco=cindex2
c_index_plco$X1=c("MetaPRS","LC-DAI-PRS","LC-GRAFF-PRS","LC-JIA-PRS","LC-SHI-PRS","LC-ZHANG-PRS",
                  "LC-FRITSCHE-14","LC-FRITSCHE-19","LC-HUNG-35","LC-HUNG-128")
fig.p2=ggplot(c_index_plco, aes(x = X1, y = cindex, fill = X1, color = "black")) + 
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
  ggtitle('PLCO')+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none")
fig.p2





#####(3)figure 4-C
rm(list = ls())
gc()

library(splines)
library(lattice)
library(rms)
library(ggplot2)
library(survminer)
library(survival)
library(data.table)

load("/Public/mzm/metaPRS/data/final.data/PLCO_baseline_20231014.Rdata")  
load('/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata')

# standardization
summary(data$metaPRS) #mean=0.157561 
summary(data_PLCO$metaPRS) #mean=7.453 
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
data_PLCO$metaPRS.sd=data_PLCO$metaPRS2/sd(data_PLCO$metaPRS2)

pcubic=datadist(data_PLCO)
options(datadist='pcubic')

fit_3=cph(Surv(lung_exitdays,lung_cancer) ~ rcs(metaPRS.sd,3)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
            chip+age+sex, data = data_PLCO,x=T,y=T)
fit_4=cph(Surv(lung_exitdays,lung_cancer) ~ rcs(metaPRS.sd,4)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
            chip+age+sex, data = data_PLCO,x=T,y=T)
fit_5=cph(Surv(lung_exitdays,lung_cancer) ~ rcs(metaPRS.sd,5)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
            chip+age+sex, data = data_PLCO,x=T,y=T)

AIC(fit_3) 
AIC(fit_4)
AIC(fit_5)

# The ref is adjusted according to y-hat =1, specifying a reference value to package the data again, OR=1
Pre_HR <-rms::Predict(fit_3,metaPRS.sd,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
ggplot(Pre_HR)
# Select the metaPRS.sd value corresponding to y-hat=1
a=subset(Pre_HR,round(Pre_HR$yhat,0)==1,select=c('yhat','metaPRS.sd')) #0.31032165

ddist <- datadist(data_PLCO)
ddist$limits$metaPRS.sd[2] <- 0.31032165
options(datadist="ddist")

fit_3.2=cph(Surv(lung_exitdays,lung_cancer) ~ rcs(metaPRS.sd,3)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
              chip+age+sex, data = data_PLCO,x=T,y=T)

anova(fit_3) ##P-value on the figure
pred2=Predict(fit_3.2,metaPRS.sd,fun=exp,ref.zero=T)

### figure
PP1=ggplot(pred2,ylim=c(0,5))+
  geom_ribbon(aes(metaPRS.sd,ymin = lower, ymax = upper),linetype='blank',linewidth=0.7,fill="lightskyblue1")+theme_bw()+
  geom_line(aes(x=metaPRS.sd ,y=yhat),linetype=1,linewidth=1,color="#0099CC")+
  geom_hline(yintercept=1, linetype=2,color="#AD002A99",linewidth=1)+
  labs(x="metaPRS (per sd increased)", y="HR of lung cancer")+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text=element_text(size=10,colour='black'),
        axis.title.x=element_text(size=12,vjust=-3),axis.title.y=element_text(size=12,vjust=5, angle = 90),plot.margin = unit(c(1, 1, 1, 1),"cm"))
PP1

# HR on the figure
output=function(x){   x <- summary(x)
HR = round(x$coef[1,2],2)
HR.confint.lower <- round(x$conf.int[1,"lower .95"],2) 
HR.confint.upper <- round(x$conf.int[1,"upper .95"],2)
p.value<-signif(x$coefficients[1,"Pr(>|z|)"],4)
p.value1=round(p.value,3)
res<-c(HR,HR.confint.lower, HR.confint.upper,p.value1)
return(res)}
fit_2=coxph(Surv(lung_exitdays,lung_cancer) ~ metaPRS.sd+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
              chip+age+sex, data = data_PLCO)
res2=output(fit_2);res2  




#####(4)figure 4-D
rm(list = ls())
gc()

library(ggplot2)
library(survminer)
library(survival)
library(data.table)

load("/Public/mzm/metaPRS/data/final.data/PLCO_baseline_20231014.Rdata") 
load('/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.RData')

data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
data_PLCO$metaPRS.sd=data_PLCO$metaPRS2/sd(data_PLCO$metaPRS2)

data$metaPRS.sd=data$metaPRS/sd(data$metaPRS)
quq=quantile(data$metaPRS.sd,seq(0.2,1,0.2))
quq

attach(data_PLCO)
data_PLCO$group=as.factor(ifelse(metaPRS.sd<=quq[1],1,ifelse(metaPRS.sd>quq[1] & metaPRS.sd<quq[4],2,3)))
detach(data_PLCO)

data_PLCO$lung_exitdays=data_PLCO$lung_exitdays/365.25

output1=function(x){   x <- summary(x)
HR = round(x$coef[c(1:2),2],2)
HR.confint.lower <- round(x$conf.int[c(1:2),"lower .95"],2) 
HR.confint.upper <- round(x$conf.int[c(1:2),"upper .95"],2)
p.value<-signif(x$coefficients[c(1:2),"Pr(>|z|)"],3)
p.value1=round(p.value,3)
res<-cbind(HR,HR.confint.lower, HR.confint.upper,p.value1)
return(res)}

fit_32=coxph(Surv(lung_exitdays,lung_cancer) ~ group+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               chip+age+sex, data = data_PLCO)
res2=output1(fit_32);res2  #HR and 95%CI in figure


####数据标准化
dat4=subset(data_PLCO,select=c('group', 'sex', 'age',
                               'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'chip'))
dat_df=with(dat4,
            data.frame(lung_exitdays=rep(10,108),
                       lung_cancer=rep(1,108),
                       sex=rep(round(median(sex,na.rm=TRUE),0),108),
                       age=rep(c(42:77),3), 
                       PC1=rep(mean(PC1,na.rm=TRUE),108),
                       PC2=rep(mean(PC2,na.rm=TRUE),108),
                       PC3=rep(mean(PC3,na.rm=TRUE),108),
                       PC4=rep(mean(PC4,na.rm=TRUE),108),
                       PC5=rep(mean(PC5,na.rm=TRUE),108),
                       PC6=rep(mean(PC6,na.rm=TRUE),108),
                       PC7=rep(mean(PC7,na.rm=TRUE),108),
                       PC8=rep(mean(PC8,na.rm=TRUE),108),
                       PC9=rep(mean(PC9,na.rm=TRUE),108),
                       PC10=rep(mean(PC10,na.rm=TRUE),108),
                       chip=rep(round(mean(chip,na.rm=TRUE),0),108),
                       group=c(rep('1',36),rep('2',36),rep('3',36))))

##计算生存时间
p2 <- predict(fit_32, newdata = dat_df, type="expected")  

#plot Figures#
plotdata=cbind(dat_df,p2)
names(plotdata)[17]=c("abosulate_rate")
plotdata$abosulate_rate=plotdata$abosulate_rate*100

p2<-ggplot(plotdata,aes(x=age,y=abosulate_rate, colour=group,linetype=group))+
  geom_line(size=0.8)+
  labs(x="Age",y="Absolute risk,%")+
  theme(panel.grid.major=element_line(colour=NA))+
  scale_linetype_manual(values=c('solid','solid','solid'))+
  scale_color_manual(values = c('#86C27D',"#4A6990FF","#DB423E"),labels=c("Low genetic risk",
                                                                          "Intermediate genetic risk: 1.52(1.36-1.70)",'High genetic risk:2.31(2.01-2.65)'))+
  theme(legend.background = element_blank(),panel.background=element_rect(fill='transparent'))+
  theme(panel.border=element_rect(color = "grey", fill = NA, linewidth = 0.8))+
  theme(legend.position = c(0.1, .95),
        legend.justification = c("left", "top"),panel.grid.major = element_line(colour = "grey", linetype = "dashed"))

p2
