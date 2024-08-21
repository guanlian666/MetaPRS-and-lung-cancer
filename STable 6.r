rm(list = ls())
gc()

load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_20231018.Rdata")  #data
names(data)

library(survival)
library(broom)
library(gtsummary)


#### logistic
#smoking initiation (SI) status
data$smoking_initiation=ifelse(data$Smoking_status==0,0,ifelse(data$Smoking_status %in% c(1,2),1,NA))
with(data,table(smoking_initiation,Smoking_status))

out_log=function(x){beta=round(coef(x)[[2]],3)
                 se=summary(x)$coefficients[2,2]
				 re=c(beta,se)  
                 beta.se=paste0(re[[1]],"(",re[[2]],")")   
                 or.ci=paste0(round(exp(re[[1]]),2),"(",round(exp(re[[1]]-1.96*re[[2]]),2),"-",round(exp(re[[1]]+1.96*re[[2]]),2),")")   
				 p=round(tidy(x)$p.value[[2]],3)
				 res=data.frame('Beta.se'=beta.se,
								'OR'=or.ci,
								'P'=p)
				 return(res)
				 }			 
				 
data$SI_PRS_persd=data$SI_PRS/sd(data$SI_PRS)
mo=glm(smoking_initiation~SI_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=binomial(),data=data)

re_smoke=out_log(mo)


#age of smoking initiation  
out_lmm=function(x){beta=round(coef(x)[[2]],3)
                 se=summary(x)$coefficients[2,2]
				 re=c(beta,se)  
                 beta.se=paste0(re[[1]],"(",re[[2]],")")  
                 or.ci=NA 
				 p=round(tidy(x)$p.value[[2]],3)
				 res=data.frame('Beta.se'=beta.se,
								'OR'=or.ci,
								'P'=p)
				 return(res)
				 }		


data$AI_PRS_persd=data$AI_PRS/sd(data$AI_PRS)				
mo2=glm(Age_started_smoking~AI_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_AI=out_lmm(mo2)


#Average cigarettes/day	
data$CPD_PRS_persd=data$CPD_PRS/sd(data$CPD_PRS)
mo2=glm(Number_of_cigarettes_daily~CPD_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_ACD=out_lmm(mo2)


# Smoking_cessation
data$SC_PRS_persd=data$SC_PRS/sd(data$SC_PRS)
mo2=glm(Smoking_cessation~SC_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=binomial(),data=data)
re_SC=out_log(mo2)



#### lung function
# FEV1
data$FEV1_PRS_persd=data$FEV1_PRS/sd(data$FEV1_PRS)
data$FEV1_max=as.numeric(data$FEV1_max)
mo2=glm(FEV1_max~FEV1_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_FEV1=out_lmm(mo2)

# FVC_PRS
data$FVC_PRS_persd=data$FVC_PRS/sd(data$FVC_PRS)
data$FVC_max=as.numeric(as.character(data$FVC_max))
mo2=glm(FVC_max~FVC_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_FVC=out_lmm(mo2)


# FEV1FVC
data$FEV1FVC_PRS_persd=data$FEV1FVC_PRS/sd(data$FEV1FVC_PRS)
data$FEV1FVC=as.numeric(data$FEV1_max)/as.numeric(data$FVC_max)
mo2=glm(FEV1FVC~FEV1FVC_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_FEV1FVC=out_lmm(mo2)


#PEF
data$PEF_PRS_persd=data$PEF_PRS/sd(data$PEF_PRS)
data$PEF_max=as.numeric(data$PEF_max)
mo2=glm(PEF_max~PEF_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_PEF=out_lmm(mo2)



#### Other body measurements
#height
data$HEIGHT_PRS_persd=data$HEIGHT_PRS/sd(data$HEIGHT_PRS)
mo2=glm(Standing_height~HEIGHT_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_height=out_lmm(mo2)


#BMI 
data$BMI_PRS_persd=data$BMI_PRS/sd(data$BMI_PRS)
mo2=glm(BMI~BMI_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_BMI=out_lmm(mo2)



#eduction
data$EDU_PRS_persd=data$EDU_PRS/sd(data$EDU_PRS)
mo2=glm(educate_2~EDU_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=gaussian(),data=data)
re_edu=out_lmm(mo2)



# Family history of lung cancer
data$LCHIST_PRS_persd=data$LCHIST_PRS/sd(data$LCHIST_PRS)
mo=glm(FH_Lung~LCHIST_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
         PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,family=binomial(),data=data)
re_LCHIST=out_log(mo)
 

	  
#### cox
output=function(x){   x <- summary(x)
 beta=round(x$coef[1,1],3)
 se=x$coef[1,3]
 beta.se=paste0(beta,"(",se,")")
 HR = round(x$coef[1,2],2)
 HR.confint.lower <- round(x$conf.int[1,"lower .95"],2) 
 HR.confint.upper <- round(x$conf.int[1,"upper .95"],2)
 hr.ci=paste0(HR,"(", HR.confint.lower ,"-",HR.confint.upper,")")
 p.value<-signif(x$coefficients[1,"Pr(>|z|)"],4)
 p.value1=round(p.value,3)
 res<-c(beta.se,hr.ci, p.value1)
 return(res)}

### AD_lung
data$LUAD_PRS_persd=data$LUAD_PRS/sd(data$LUAD_PRS)
x=coxph(Surv(difftime,AD)~LUAD_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data)
re.AD=output(x)                 


### SC_lung
data$LUSC_PRS_persd=data$LUSC_PRS/sd(data$LUSC_PRS)
x=coxph(Surv(difftime,event=SC)~LUSC_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data)
re.sc=output(x)     

		
### SLC_lung
data$SCLC_PRS_persd=data$SCLC_PRS/sd(data$SCLC_PRS)
x=coxph(Surv(difftime,event=SLC)~SCLC_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data)
re.slc=output(x)      



###asthma
load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_asthma.rdata") 
data_asthma$ASTHMA_PRS_persd=data_asthma$ASTHMA_PRS/sd(data_asthma$ASTHMA_PRS)
x=coxph(Surv(asthma_time,asthma_ICD10)~ASTHMA_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data_asthma)
re.asthma=output(x) 


###COPD  
load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_COPD.rdata")
data_COPD$COPD_PRS_persd=data_COPD$COPD_PRS/sd(data_COPD$COPD_PRS)
x=coxph(Surv(COPD_time,COPD_ICD10)~COPD_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data_COPD)
re.COPD=output(x)                


###IPF  
load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_IPF.rdata")
data_IPF$IPF_PRS_persd=data_IPF$IPF_PRS/sd(data_IPF$IPF_PRS)
x=coxph(Surv(IPF_time,IPF_ICD10)~IPF_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data_IPF)
re.IPF=output(x) 
    
   
###ILD  
load("/Public/mzm/metaPRS/data/final.data/UKB_baseline_ILD.rdata")
data_ILD$ILD_PRS_persd=data_ILD$ILD_PRS/sd(data_ILD$ILD_PRS)
x=coxph(Surv(ILD_time,ILD_ICD10)~ILD_PRS_persd+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+
        PCA7+PCA8+PCA9+PCA10+chip+Age_when_attended_assessment_centre+Sex,data=data_ILD)
re.ILD=output(x)


res=rbind(re.AD,re.sc,re.slc,re_smoke,re_AI,re_ACD,re_SC,re_FEV1,re_FVC,re_FEV1FVC,re_PEF,
          re.COPD,re.IPF,re.ILD,re.asthma,re_height,re_BMI,re_edu,re_LCHIST)

name=c('LUAD','LUSQ','SCLC','SI','AI','CPD','SC','FEV1','FVC','FEV1/FVC','PEF',
       'COPD','IPF','ILD','asthma','height','BMI','Education','Family history of lung cancer')

re2=cbind(name,res)	   
write.csv(re2,file='STable6.csv')

