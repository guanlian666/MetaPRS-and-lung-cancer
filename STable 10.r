library(interactionR)
library(epiDisplay)
load("PLCO_baseline_20231014.Rdata")  

#######################################################################################################################
da0<- subset(data_PLCO,(educat1!=9) & (COPD!=9) & (lung_fh!=9) & (race7!=99) & (ph_any_cancer!=9))  

da0$race7<-as.factor(da0$race7)
da0$educat1<-as.factor(da0$educat1)
da0$COPD<-as.factor(da0$COPD)
da0$ph_any_cancer<-as.factor(da0$ph_any_cancer)
da0$fh_cancer<-as.factor(da0$fh_cancer)
da0$cig_stat<-as.factor(da0$cig_stat)

da0[da0$is.smoke==0,]$cig_stop1<-10


#Smoking_initiation
da0$Smoking_initiation=ifelse(da0$is.smoke==1 & da0$Smoking_cessation==1,2,
                 ifelse(da0$is.smoke==1 & da0$Smoking_cessation==0,1,
				 ifelse(da0$is.smoke==0,0,NA)))
da0$Smoking_initiation=as.factor(da0$Smoking_initiation)


####interaction
#age				 
model<-glm(lung_cancer~age*metaPRS+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
summary(model)


##Ethnicity
model<-glm(lung_cancer~age+metaPRS+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
model2<-glm(lung_cancer~age+metaPRS*race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
lrtest(model2,model) 


# Education
model2<-glm(lung_cancer~age+metaPRS*educat1+race7+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
lrtest(model2,model)


# BMI		 
model<-glm(lung_cancer~age+race7+educat1+BMI*metaPRS+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
summary(model)



# COPD
model<-glm(lung_cancer~age+metaPRS+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD*metaPRS+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
lrtest(model2,model)


# Personal history of cancer
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD+ph_any_cancer*metaPRS+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
lrtest(model2,model)



# Family history of lung cancer
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh*metaPRS+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
lrtest(model2,model)

 
# Smoking initiation
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation*metaPRS+cigpd_f1+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
lrtest(model2,model)


# Cigarettes per day
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1*metaPRS+cig_years1+cig_stop1,family = binomial(link="logit"),data=da0)
summary(model2)
   

# Smoking duration
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1*metaPRS+cig_stop1,family = binomial(link="logit"),data=da0)
summary(model2)



# Smoking quit time
model2<-glm(lung_cancer~age+race7+educat1+BMI+
             COPD+ph_any_cancer+
             lung_fh+Smoking_initiation+cigpd_f1+cig_years1+cig_stop1*metaPRS,family = binomial(link="logit"),data=da0)
summary(model2)




