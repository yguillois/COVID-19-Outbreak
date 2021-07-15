
##### PACKAGES #####

library(tidyverse)
library(MASS)

# overdispersion
library(AER)

# Confint
library(emmeans)

# Effects
library(effects)
library(ggeffects)

##### LOADING #####
MP <- read_delim("S:/DIRE/CIRE_BRE/Cluster Abattoir/KERMENE/Dossier pour analyse rétrospective du cluster/MP.csv",
";", escape_double = FALSE, trim_ws = TRUE)

##### ANALYSES #####

#### VARIABLE CREATION ####

## Age ##
MP$Age<-2020 - MP$year_of_birth
MP$Age10<-(2020 - MP$year_of_birth)/10
MP<-MP %>%
mutate(age_group = cut(Age, breaks = c(0,18,30,40,50,60,150), right = FALSE))

## Sex ##
MP$Sex<-forcats::fct_relevel(MP$Sex,"Women")

# Country
MP$Country<-forcats::fct_relevel(MP$Country,"FRANCAISE")

# France vs other
MP$France<-ifelse(MP$Country=="FRANCAISE","FRANCAISE","ETRANGERE")
MP$France<-forcats::fct_relevel(MP$France,"FRANCAISE")

## Language ##

# French language
MP$Fr_language<-ifelse(MP$language=="FRANCOPHONE","FRANCOPHONE","OTHER")
MP$Fr_language<-forcats::fct_relevel(MP$Fr_language,"FRANCOPHONE")

## Employer##
MP<-MP%>%

mutate(Employer2=if_else(Employer %in% c("Regular","Veterinary"),"Vet-Reg",Employer))%>%         
mutate(Employer2 = forcats::fct_relevel(Employer2,"Vet-Reg","Temporary","Subcontractor"))%>%

mutate(Employer_cor2=if_else(Employer_cor %in% c("Regular","Veterinary"),"Vet-Reg",Employer_cor))%>%         
mutate(Employer_cor2 = forcats::fct_relevel(Employer_cor2,"Vet-Reg","Temporary","Subcontractor"))%>%

## Work hours##
mutate(Shift2=if_else(Shift %in% c("Night","Flextime","In the day"),"Ref",Shift))%>%
mutate(Shift2 = forcats::fct_relevel(Shift2,"Ref","Morning","Afternoon"))%>%

## Activities ##      
mutate(Activity2=if_else(Activity %in% c("S","T","Other"),"Ref",Activity))%>%
mutate(Activity2=forcats::fct_relevel(Activity2,"Ref","C"))%>%
   

## Campaign (at least, one campaign)
mutate(Campaign=case_when(Campaign_15 %in% c("Yes") | Campaign_19 %in% c("Yes") | Campaign_25 %in% c("Yes") | Campaign_26 %in% c("Yes")~1,TRUE~0))%>%


## Posititive RT_PCR during one campaign ##
mutate(Pos_Campaign = case_when (Res15 %in% c("POSITIF")~"POS", 
                             Res19 %in% c("POSITIF")~"POS",
                             Res25 %in% c("POSITIF")~"POS",
                             Res26 %in% c("POSITIF")~"POS",           
                             TRUE~"NOT POS"))%>%

## car pooling or sharing accomodation ##
mutate(cp_sa = case_when(carpooling %in% c("Yes")~"Yes",
                         accomodation %in% c("Yes")~"Yes",
                         Forms %in% c("Yes")~"No",
                         TRUE ~ NA_character_))

#### DESCRIPTION OF THE STUDY POPULATION ####

## Sex and Age ##
addmargins(table(MP$Sex,useNA="always"))
prop.table(table(MP$Sex,useNA="always"))
prop.table(table(MP$Sex))
summary(MP$Age)

## Employers ##
addmargins(table(MP$Employer,useNA="always"))
addmargins(table(MP$Employer))
prop.table(table(MP$Employer))

table(MP$Employer_det[MP$Employer=="Subcontractor"])
table(MP$Employer_det[MP$Employer=="Subcontractor"],MP$Activity[MP$Employer=="Subcontractor"])
addmargins(prop.table((table(MP$Employer_det[MP$Employer=="Subcontractor"],MP$Activity[MP$Employer=="Subcontractor"]))))

## Country of birth ##
addmargins(table(MP$France))
prop.table(table(MP$France))
addmargins(table(MP$EEW))
prop.table(table(MP$EEW))

## Tested workers ##
addmargins(table(MP$T))
prop.table(table(MP$T))

with(MP,by(Age,T,summary))
sr<-addmargins(table(MP$Sex,MP$T))
Sex_ratio_MF_non_tested<-(sr[2,1]/sr[1,1])
Sex_ratio_MF_tested<-(sr[2,2]/sr[1,2])
print(paste("Sex ratio tested :",Sex_ratio_MF_tested,sep=" " ))
print(paste("Sex ratio non tested :",Sex_ratio_MF_non_tested,sep=" " ))

## Tested workers and Employers ##
addmargins(table(MP$Employer,MP$T))
addmargins(prop.table(table(MP$Employer,MP$T),margin=1),margin=2)

## Sampled during one campaign ##
addmargins(table(MP$Campaign[MP$T==1]))
prop.table(table(MP$Campaign[MP$T==1]))

## Workers tested repeatedly (including hospital and outpatient sampling) ##
tr<-table(MP$Num_tests)
print(tr)
print(paste("%=",100*(tr[3]+tr[4])/sum(MP$T),sep=" "))

addmargins(table(MP$A_workshop[MP$Num_tests>1]))
prop.table(table(MP$A_workshop[MP$Num_tests>1]))

#### DESCRIPTION OF OCCUPATIONAL CASES ####

## Tests and campaigns ##
addmargins(table(MP$case))
addmargins(table(MP$Pos_Campaign [MP$case==1]))
addmargins(prop.table(table(MP$Pos_Campaign [MP$case==1])))

addmargins(table(MP$date [MP$case==1]))
addmargins(table(MP$date,MP$case,useNA="always"),margin=1)

## Age ##
cases<-MP%>%
   filter(case==1)

summary(cases$Age)

## Sex_ratio ##
sr_case<-table(cases$Sex)
print(paste("sex-ratio=",sr_case[2]/sr_case[1], sep=" "))

## Distance ##
summary(cases$D)

## Symptoms ##
addmargins(table(cases$Symptoms,useNA="always"))
prop.table(table(cases$Symptoms))

## Hospitalized and ICU ##
# Numbers
table(cases$hosp)
prop.table(table(cases$hosp))
summary(cases$Age[cases$hosp==1])
table(cases$ICU)
# Activity
table(cases$Activity[cases$hosp==1])
table(cases$d1hosp)

## Employer ##
addmargins(table(cases$Employer))
prop.table(table(cases$Employer))

## Shift ##
addmargins(table(cases$Shift,useNA="always"))
prop.table(table(cases$Shift))

## Attack rates and distributions in workshops##

# attack rate in the study population
table(MP$case)
prop.table(table(MP$case))

# attack rate in the cutting department
addmargins(table(MP$case,MP$Activity,useNA="always"))
prop.table(table(MP$case,MP$Activity,useNA="always"),margin=2)

# attack rates in A_workshop and primary cutting workshop 
addmargins(table(MP$case,MP$A_workshop,useNA="always"))
prop.table(table(MP$case,MP$A_workshop,useNA="always"),margin=2)

# distributions 
addmargins(table(cases$Activity))
prop.table(table(cases$Activity))
addmargins(table(cases$A_workshop))
prop.table(table(cases$A_workshop))

## Positivity rate in A_workshop ##
AW<-MP%>%
   filter(A_workshop=="A_workshop")%>%
   dplyr::select(ID,case,dout_hos,Res15,Res19,Res26,Campaign_15,Campaign_19,Campaign_26,A_workshop)

AWc<-AW%>%
   filter(case==1)

# 1st screening campaign
AWc15<-AW%>%
   filter(Res15=="POSITIF" & case==1)

posrat1<-table(AW$Res15[AW$Campaign_15=="Yes"])
print(posrat1)
print(paste("positivity rate =",(posrat1[4]/(posrat1[4]+posrat1[2]))))



# 19 may campaign campaign
table(MP$Res19[MP$Campaign_19=="Yes" & MP$Campaign_15=="No" & MP$A_workshop=="A_workshop" ],useNA="always")

AWc19<-AW%>%
   filter(Res15 !="POSITIF" | Campaign_15=="No")%>%
   filter(Res19=="POSITIF")


# Cases detected during the 2 first campaigns in the A-Workshop
AWc15_19<-bind_rows(AWc15,AWc19)%>%
   dplyr::select(ID)%>%
   mutate(ind=1)


# Supplementary cases detected during the 26 May campaign
AWc26<-left_join(AW,AWc15_19,by="ID") %>%
   filter(Campaign_26=="Yes")%>%
   filter(is.na(ind))

table(AWc26$Res26,useNA="always")
prop.table(table(AWc26$Res26,useNA="always"))

AWc26<-AWc26 %>%
   filter(Res26=="POSITIF")%>%
   dplyr::select(ID)%>%
   mutate(ind=1)

AWcases<-bind_rows(AWc15_19,AWc26)
   
# Outpatient sampling
AWcases_2<-left_join(AWc,AWcases)%>%
   filter(is.na(ind))

table(AWcases_2$dout_hos,useNA = "always")

## Foreign-born cases ##
addmargins(table(cases$France))
prop.table(table(cases$France))

table(MP$France, MP$case)
prop.table(table(MP$France, MP$case),margin=2)

table(cases$France, cases$EEW, useNA="always")
prop.table(table(cases$EEW[cases$France=="ETRANGERE"]))

table(cases$France, cases$Fr_language, useNA="always")
prop.table(table(cases$Fr_language[cases$France=="ETRANGERE"]))

table(cases$Activity [cases$EEW=="EEW"], useNA="always")
prop.table(table(cases$Activity [cases$EEW=="EEW"], useNA="always"))

table(cases$Employer [cases$EEW=="EEW"], useNA="always")
prop.table(table(cases$Employer [cases$EEW=="EEW"], useNA="always"))

## Contact-tracing forms ##
table(cases$Forms,useNA="always")
prop.table(table(cases$Forms,useNA="always"))
# foreign-born
table(cases$France[cases$Forms=="Yes"],useNA="always")
# car-pooling
addmargins(table(cases$carpooling [cases$Forms=="Yes"], cases$Forms[cases$Forms=="Yes"],useNA="always"),margin=1)
prop.table(table(cases$carpooling [cases$Forms=="Yes"], cases$Forms[cases$Forms=="Yes"],useNA="always"),margin=2)
# sharing accomodation
addmargins(table(cases$accomodation [cases$Forms=="Yes"], cases$Forms[cases$Forms=="Yes"],useNA="always"),margin=1)
prop.table(table(cases$accomodation [cases$Forms=="Yes"], cases$Forms[cases$Forms=="Yes"],useNA="always"),margin=2)

# car pooling or sharing accomodation
addmargins(table(cases$cp_sa, cases$EEW))
addmargins(prop.table(table(cases$cp_sa, cases$EEW),margin=2),margin=1)

chisq.test(cases$cp_sa,cases$EEW) 
# chisq.test(cases$cp_sa,cases$EEW,simulate.p.value = TRUE)

#### UNIVARIATE ANALYSES ####

## Sex ##
addmargins(table(MP$Sex[MP$T==1]))
prop.table(table(MP$Sex [MP$T==1]))
addmargins(table(MP$Sex,MP$case))
prop.table(table(MP$Sex,MP$case),margin=2)


model <- glm(case ~ Sex, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP Sex et son IC
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))

# p de Sex
smodel<-summary(model)
print(smodel)

## Age ##
addmargins(table(MP$age_group[MP$T==1]))
prop.table(table(MP$age_group [MP$T==1]))
addmargins(table(MP$age_group,MP$case))
prop.table(table(MP$age_group,MP$case),margin=2)

model <- glm(case ~ Age10, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP Age10 et son IC
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))

# p de Age10
smodel<-summary(model)
print(smodel)

## Place of birth - Foreign-born vs born in France ##
addmargins(table(MP$France[MP$T==1]))
prop.table(table(MP$France [MP$T==1]))
addmargins(table(MP$France,MP$case))
prop.table(table(MP$France,MP$case),margin=2)

model <- glm(case ~ France, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP France-born et son IC
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))

# p de France-born
smodel<-summary(model)
print(smodel)

## Place of birth - EEW ##
addmargins(table(MP$EEW[MP$T==1]))
prop.table(table(MP$EEW[MP$T==1]))
addmargins(table(MP$EEW,MP$case))
prop.table(table(MP$EEW,MP$case),margin=2)

model <- glm(case ~ EEW, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP EEW et son IC
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))

# p de EEW
smodel<-summary(model)
print(smodel)

## Employer ##
addmargins(table(MP$Employer[MP$T==1]))
prop.table(table(MP$Employer[MP$T==1]))
addmargins(table(MP$Employer,MP$case))
prop.table(table(MP$Employer,MP$case),margin=2)

addmargins(table(MP$Employer2[MP$T==1]))

model <- glm(case ~ Employer2, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP Employer2 et son IC
# RP temporary workers
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))
# RP subcontractors
print(c(exp(model$coefficients[3]),exp(u[3,1]),exp(u[3,2])))

# p de Employer2
smodel<-summary(model)
print(smodel)

## Activity and A_workshop ##
addmargins(table(MP$Activity[MP$T==1]))
prop.table(table(MP$Activity[MP$T==1]))
addmargins(table(MP$Activity,MP$case))
prop.table(table(MP$Activity,MP$case),margin=2)

addmargins(table(MP$A_workshop[MP$T==1]))
prop.table(table(MP$A_workshop[MP$T==1]))
addmargins(table(MP$A_workshop,MP$case))
prop.table(table(MP$A_workshop,MP$case),margin=2)

addmargins(table(MP$Activity2[MP$T==1]))

model <- glm(case ~ Activity2, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP Activity2 et son IC
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))

# p de Activity2
smodel<-summary(model)
print(smodel)


## Shift ##
addmargins(table(MP$Shift[MP$T==1]))
prop.table(table(MP$Shift[MP$T==1]))
addmargins(table(MP$Shift,MP$case))
prop.table(table(MP$Shift,MP$case),margin=2)

addmargins(table(MP$Shift2[MP$T==1]))


model <- glm(case ~ Shift2, family=poisson(link=log), data=MP)
print(model)

u<-confint.default(model,level=0.95)
# RP Shift2 et son IC
# Morning
print(c(exp(model$coefficients[2]),exp(u[2,1]),exp(u[2,2])))
# Afternoon
print(c(exp(model$coefficients[3]),exp(u[3,1]),exp(u[3,2])))

# p de Shift2
smodel<-summary(model)
print(smodel)

# timetable and employer2 
addmargins(table(MP$Shift2[MP$T==1],MP$Employer2[MP$T==1],useNA="always"))
tt<-addmargins(prop.table(table(MP$Shift2[MP$T==1],MP$Employer2[MP$T==1],useNA="always"),margin=2),margin=1)
print(paste("Vet-Reg",1-tt[4,1],"Temporary",1-tt[4,2],"Subcontractors",1-tt[4,3], sep=" "))
#### MULTIVARIABLE MODEL ####
### Selection of the variables ###
MP2<-MP%>%
   dplyr::select(case,Sex,Age,France,EEW,Employer2,Activity2,Shift2)%>%
   dplyr::filter(!is.na(case) & !is.na(Sex) & ! is.na(Activity2) &!is.na(Age) & !is.na(Shift2) & !is.na(France) & !is.na(EEW) & !is.na(Employer2) )%>%
   mutate(Age=as.factor(Age),case=as.vector(case))



model_c<-glm(case ~ Sex + Activity2+ Age + Shift2+ France + EEW+ Employer2 , family=poisson(link=log), data=MP2, na.action = na.omit)
print(model_c)
   
model_comp2 <-stepAIC(model_c,direction ="backward", k =log(nrow(MP2)))
summary(model_comp2)

### Interaction terms ###
MP3<-MP%>%
   dplyr::select(case,EEW,Employer2,Activity2)%>%
   dplyr::filter(!is.na(case) & ! is.na(Activity2) & !is.na(EEW) & !is.na(Employer2) )%>%
   mutate(case=as.vector(case))

model_2c <- glm(case ~ (EEW + Employer2 + Activity2)^2, family=poisson(link=log), data=MP3, na.action = na.omit )

summary(model_2c)
summary(model_2c, correlation=T)

### selected model ###
sel_model<-glm(case ~ Employer2 + EEW*Activity2, family=poisson(link=log), data=MP3,na.action = na.omit )
   
# Ajustement 
print(pchisq(sel_model$deviance,sel_model$df.residual,lower.tail =FALSE))
khi2<-sum((MP3$case-sel_model$fitted.values)^2/sel_model$fitted.values)
print(pchisq(khi2,sel_model$df.residual,lower.tail =FALSE))

# Overdispersion . 
ratio_overdispersion <-(khi2/sel_model$df.residual)
print(ratio_overdispersion)

dispersiontest(sel_model)

um_sel_model<-confint(sel_model,level=0.95)

# RR temporary
print(c(exp(sel_model$coefficients[2]),exp(um_sel_model[2,1]),exp(um_sel_model[2,2])))
# RR subcontractor
print(c(exp(sel_model$coefficients[3]),exp(um_sel_model[3,1]),exp(um_sel_model[3,2])))
# EEW, outside the cutting department
print(c(exp(sel_model$coefficients[4]),exp(um_sel_model[4,1]),exp(um_sel_model[4,2])))
# Other workers, Cutting department
print(c(exp(sel_model$coefficients[5]),exp(um_sel_model[5,1]),exp(um_sel_model[5,2])))
# EEW, Cutting department
cov<-summary(sel_model)$cov.scaled
varEEW<-cov[4,4]
varEEWcut<-cov[6,6]
covar<-cov[4,6]
var<-varEEW+varEEWcut+2*covar
rr<-exp(sel_model$coefficients[4]+sel_model$coefficients[6])
rrinf<-exp(sel_model$coefficients[4]+sel_model$coefficients[6]-2*(var^.5))
rrsup<-exp(sel_model$coefficients[4]+sel_model$coefficients[6]+2*(var^.5))
print(c(rr,rrinf,rrsup))

predictorEffects(sel_model)
ris<-as.data.frame(predictorEffects(sel_model))

# Vet_Reg
print(ris$Employer2[1,1])
print(c(100*ris$Employer2[1,2],100*ris$Employer2[1,4],100*ris$Employer2[1,5]))
# Temporary
print(ris$Employer2[2,1])
print(c(100*ris$Employer2[2,2],100*ris$Employer2[2,4],100*ris$Employer2[2,5]))
# Subcontractors
print(ris$Employer2[3,1])
print(c(100*ris$Employer2[3,2],100*ris$Employer2[3,4],100*ris$Employer2[3,5]))
# EEW in the cutting department
print(ris$EEW[1,1])
print(c(100*ris$EEW[4,3],100*ris$EEW[4,5],100*ris$EEW[4,6]))
# EEW outside the cutting department
print(c(100*ris$EEW[2,3],100*ris$EEW[2,5],100*ris$EEW[2,6]))
# Other workers in the cutting department
print(c(100*ris$EEW[3,3],100*ris$EEW[3,5],100*ris$EEW[3,6]))
# Other workers outside the cutting department
print(c(100*ris$EEW[1,3],100*ris$EEW[1,5],100*ris$EEW[1,6]))