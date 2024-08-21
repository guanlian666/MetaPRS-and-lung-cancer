rm(list=ls())
gc()

setwd('/Public/mzm/metaPRS/data/final.data')
load('PLCO_baseline_20231014.Rdata')   #data_PLCO

names(data_PLCO)
data_PLCO$cig_stop2=ifelse(is.na(data_PLCO$cig_stop1)==T&data_PLCO$cig_stat==0,10,data_PLCO$cig_stop1)
data_PLCO$arm1=ifelse(data_PLCO$arm==2,0,data_PLCO$arm)   #0 control, 1 intervention


###metaPRS group: Divided into 3 layers according to the top 20% and bottom 80% of metaPRS in UKB cohort
load('UKB_baseline_20231018.Rdata')
quq=quantile(data$metaPRS,seq(0.2,1,0.2))
quq
data_PLCO$metaPRS2=data_PLCO$metaPRS-(7.453-0.157561)
data_PLCO$PRS3=as.numeric(ifelse(data_PLCO$metaPRS2<=quq[1],0,ifelse(data_PLCO$metaPRS2>quq[1] & data_PLCO$metaPRS2<quq[4],1,2)))
table(data_PLCO$PRS3)


#### ---------------------(1)Incidence of lung cancer
LC_ins=subset(data_PLCO,data_PLCO$lung_cancer==1)   
summary(LC_ins$lung_exitage)  

# Groups by age of onset
LC_ins$ins_age_group=ifelse(LC_ins$lung_exitage>=55&LC_ins$lung_exitage<60,1,
                        ifelse(LC_ins$lung_exitage>=60&LC_ins$lung_exitage<65,1,
                               ifelse(LC_ins$lung_exitage>=65&LC_ins$lung_exitage<70,2,
                                      ifelse(LC_ins$lung_exitage>=70&LC_ins$lung_exitage<75,3,
                                             ifelse(LC_ins$lung_exitage>=75&LC_ins$lung_exitage<80,4,
                                                    ifelse(LC_ins$lung_exitage>=80,5,6))))))
# library(tidyverse)
library(dplyr)
LC_ins1=LC_ins[,c('#FID','PRS3','ins_age_group')]
dat = count(LC_ins1,ins_age_group,PRS3)
dat = dat %>% group_by(ins_age_group) %>% 
  reframe(PRS = PRS3,n = n/sum(n))
dat$PRS = factor(dat$PRS, levels = c(0,1,2),labels = c("Low","Intermediate","High"))
dat$ins_age_group <- factor(dat$ins_age_group,
                               levels = c(1,2,3,4,5),
                               labels = c("55~65", "65~70","70~75","75~80", "≥80"))
dat$n=round(dat$n,3)
head(dat)


#plot figure
library(ggplot2)
p1 <- ggplot(dat,aes(x=ins_age_group, y=n),position="stack") +
  scale_x_discrete(limits=c("55~65", "65~70","70~75","75~80", "≥80"))+
  geom_bar(aes(x=ins_age_group, y=n, fill = PRS), stat = "identity",color="black",linewidth=0.4,
           position = position_stack(reverse=TRUE), width = 0.6,data=dat)+ 
  geom_text(aes(x=ins_age_group, y=n,label = scales::percent(n)),
            color = "black",size = 3,
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c("#56B4E9", '#8DD3C7', '#FB8072'))+
  xlab('Age Group of Onset') + ylab("Distribution of Sample Proportion") +
  theme_bw()+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=15,face="plain",color="black"),
    legend.title=element_text(size=15,face="plain",color="black"),
    legend.position = "top",
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", linewidth = 0.4))+
  theme(text=element_text(family="A",size=15))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))
p1





####-------------- (2) Lung cancer deaths
LC_mor=subset(data_PLCO,data_PLCO$f_dthl==1)  
summary(LC_mor$mortality_exitage)  

# Group by age of death
LC_mor$mor_age_group=ifelse(LC_mor$mortality_exitage>=55&LC_mor$mortality_exitage<60,1,
                             ifelse(LC_mor$mortality_exitage>=60&LC_mor$mortality_exitage<65,1,
                                  ifelse(LC_mor$mortality_exitage>=65&LC_mor$mortality_exitage<70,2,
                                  ifelse(LC_mor$mortality_exitage>=70&LC_mor$mortality_exitage<75,3,
                                  ifelse(LC_mor$mortality_exitage>=75&LC_mor$mortality_exitage<80,4,
                                  ifelse(LC_mor$mortality_exitage>=80,5,6))))))
table(LC_mor$mor_age_group,useNA = 'always')


library(dplyr)
LC_mor1=LC_mor[,c('#FID','PRS3','mor_age_group')]
dat1 = count(LC_mor1,mor_age_group,PRS3)
dat1 = dat1 %>% group_by(mor_age_group) %>% 
  reframe(PRS = PRS3,n = n/sum(n))
dat1$PRS = factor(dat1$PRS, levels = c(0,1,2),labels = c("Low","Intermediate","High"))
dat1$mor_age_group <- factor(dat1$mor_age_group,
                            levels = c(1,2,3,4,5),
                            labels = c("55~65", "65~70","70~75","75~80", "≥80"))
dat1$n=round(dat1$n,3)
head(dat1)


p2 <- ggplot(dat1,aes(x=mor_age_group, y=n),position="stack") +
  scale_x_discrete(limits=c("55~65", "65~70","70~75","75~80", "≥80"))+
  geom_bar(aes(x=mor_age_group, y=n, fill = PRS), stat = "identity",color="black",linewidth=0.4,
           position = position_stack(reverse=TRUE), width = 0.6,data=dat1)+
  geom_text(aes(x=mor_age_group, y=n,label = scales::percent(n)),
            color = "black",size = 3,
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c("#56B4E9", '#8DD3C7', '#FB8072'))+
  xlab('Age Group of Death') + ylab("Distribution of Sample Proportion") +
  theme_bw()+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=15,face="plain",color="black"),
    legend.title=element_text(size=15,face="plain",color="black"),
    legend.position = "top",
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", linewidth = 0.4))+
  theme(text=element_text(family="A",size=15))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))
p2


library(patchwork)
fig=p1+p2+plot_annotation('A')
ggsave(fig,file='/Public/mzm/metaPRS/data/result/Sfig3.pdf')


