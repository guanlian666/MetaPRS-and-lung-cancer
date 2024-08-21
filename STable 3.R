R
rm(list=ls())
gc()
setwd('/Public/mzm/metaPRS/data/final.data')
################## UKB
library(data.table)
library(tableone)
library(knitr)

load('UKB_baseline_20231018.RData')
data$smoking_initiation=ifelse(data$Smoking_status==0,0,ifelse(data$Smoking_status %in% c(1,2),1,99))

a <- CreateTableOne(vars=c("Age_when_attended_assessment_centre","Sex","Ethnic_background2","BMI","educate_2","smoking_initiation","Number_of_cigarettes_daily",
                    "smoke_year","Smoking_cessation","Smoking_quit_time","family_cancer"),     
                    #Vector of variables to summarize
                    data = data,
                    includeNA = T,
                    addOverall = TRUE,
                    strata="outcome_LC",
                    factorVars=c("Ethnic_background2","educate_2","family_cancer","Sex","smoking_initiation","Smoking_cessation"))

a_csv<- print(a, 
              exact =c("Ethnic_background2","educate_2","family_cancer","Sex","smoking_initiation","Smoking_cessation"),
              smd=T, 
              showAllLevels = TRUE,
              quote = FALSE, 
              noSpaces = TRUE, 
              printToggle = FALSE)


kable(a_csv,  
      align = 'c', 
      caption = 'STable 3: Comparison of lung cancer')
write.csv(a_csv,'/Public/mzm/metaPRS/data/result/Stable3.UKB_baseline.csv')





####PLCO
load('PLCO_baseline_20231014.Rdata')   #data_PLCO
names(data_PLCO)


#旧data_PLCO$cig_stat
#0="Never Smoked Cigarettes" 
#1="Current Cigarette Smoker" 
#2="Former Cigarette Smoker" 

## new data_PLCO$cig_stat1
#0 No smoking
#1 Smoking
data_PLCO$cig_stat1=ifelse(data_PLCO$cig_stat==1|data_PLCO$cig_stat==2,1,data_PLCO$cig_stat)
table(data_PLCO$cig_stat1,useNA = 'always')


data_PLCO$Smoking_cessation=ifelse(data_PLCO$cig_stat==2,1,
                                   ifelse(data_PLCO$cig_stat==1,0,9))
table(data_PLCO$Smoking_cessation,useNA = 'always')

data_PLCO$educat2=ifelse(data_PLCO$educat1==0,9,data_PLCO$educat1)
b<- CreateTableOne(vars=c("age","sex","race7","bmi_curr","educat2","cig_stat1","cigpd_f1",
    "cig_years1", "Smoking_cessation","cig_stop1","fh_cancer","ph_any_cancer"), 
                  data = data_PLCO,
                  includeNA = T,
                  addOverall = TRUE,
                  strata="lung_cancer",
                   factorVars=c("sex","race7","educat2","cig_stat1","Smoking_cessation","fh_cancer","ph_any_cancer")) 

b_csv<- print(b, 
              exact =c("sex","race7","educat2","cig_stat1","Smoking_cessation","fh_cancer","ph_any_cancer"),
              smd=T, 
              showAllLevels = TRUE,
              quote = FALSE, 
              noSpaces = TRUE, 
              printToggle = FALSE)

kable(b_csv,  
      align = 'c', 
      caption = 'STable 3: Comparison of lung cancer')

write.csv(b_csv,'/Public/mzm/metaPRS/data/result/Stable3.PLCO_baseline.csv')







