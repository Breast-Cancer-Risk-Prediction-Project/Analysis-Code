#Table1. Descriptive Data for calibration cohort
#This program will create Table 1 for the sample used for model calibration
#This is limited to those >=20y and <=75y and excludes those missing time of BC dx

install.packages("table1")
library(table1)
#Read in data - this was created in Step1.DataCleaning.R
#Specify filepath

valdata<-read.csv("FILEPATH/bcrpp_valdata.R")

#Set up data for table
valdata$hrt_use<-factor(valdata$hrt_use, levels=c(1,2,3,4), labels=c("Never use", "Current use", "Former use", "Ever use, unk if current"))
valdata$pmh_type_any<-factor(valdata$pmh_type_any, levels=c(1,2,3,4,5,6,7), labels=c("Never", "Past, unk type", "Past, E only", "Past, Combined", "Current, unk type", "Current, E only", "Current, combined"))
valdata$alcohol_status<-factor(valdata$alcohol_status, levels=c(1,2,3,4), labels=c("Current drinker", "Former drinker", "Ever drinker (unk if current)", "Never drinker"))
valdata$alcgm_cat4<-factor(valdata$alcgm_cat4, levels=c(1,2,3,4), labels=c("None"=1, "<11g", "11-<22g", "22-<50g"))
valdata$smoking_status<-factor(valdata$smoking_status, levels=c(1,2,3,4),  labels=c("Current smoker", "Former smoker", "Ever smoker (unk if current)", "Never smoker"))
valdata$ocuse<-factor(valdata$ocuse, levels=c(1,2,3,4), labels=c("Never","Ever, unk. past or current", "Past", "Current"))
valdata$bmi_cat4<-factor(valdata$bmi_cat4, levels=c(1,2,3,4), labels=c("<18.5", "18.5-<25", "25-<30","30-50"))
valdata$parous<-factor(valdata$parous, levels=c(0,1), labels=c("Non-parous","Parous"))
valdata$case_5y<-factor(valdata$case_5y, levels=c(0,1),labels=c("No","Yes"))
valdata$case_5yinv<-factor(valdata$case_5yinv, levels=c(0,1), labels=c("No","Yes"))
valdata$case_10y<-factor(valdata$case_10y, levels=c(0,1),labels=c("No","Yes"))
valdata$case_10yinv<-factor(valdata$case_10yinv, levels=c(0,1),labels=c("No","Yes"))
valdata$case<-factor(valdata$case, levels=c(0,1),labels=c("No","Yes"))
valdata$bca_invasive<-factor(valdata$bca_invasive, levels=c(0,1),labels=c("No","Yes"))
valdata$meno_status2<-factor(valdata$meno_status2, levels=c(1,2), labels=c("Pre","Post"))
valdata$race<-factor(valdata$race, levels=c(1,2,3,4,5,6), labels=c("White", "Black", "Asian","Native Hawaiian/Pacific Islander", "American Indian, Alaska Native", "Other/multiracial"))
#Demographics
label(valdata$age)<-"Age (y)"
label(valdata$birth_year)<-"Birth year"
label(valdata$race)<-"Race"
#Anthropometrics
label(valdata$height)<-"Height (cm)"
label(valdata$bmi)<-"BMI (kg/m2)"
label(valdata$bmi_cat4)<-"BMI category"

#Reproductive related
label(valdata$agemenarche)<-"Age at menarche (y)"
label(valdata$parous)<-"Parous"
label(valdata$parity_new)<-"Parity"
label(valdata$ageflb)<-"Age at first birth (y), if parous"
label(valdata$ocuse)<-"Oral contraceptive use"
label(valdata$meno_status2)<-"Menopausal status"
label(valdata$hrt_use)<-"Hormone replacement therapy use"
label(valdata$pmh_type_any)<-"Type of HRT"

#Lifestyle vars
label(valdata$alcohol_status)<-"Alcohol use, status"
label(valdata$alcgm)<-"Alcohol (g/day)"
label(valdata$smoking_status)<-"Smoking status"
label(valdata$smoking_amt)<-"Smoking, cigs/day"

#Health history vars
label(valdata$famhx_first)<-"First-degree family history BC"
label(valdata$ever_biopsy)<-"Ever had a biopsy"
label(valdata$biopsy_num)<-"Number of biopsies"
label(valdata$bbd_history)<-"History of BBD"
label(valdata$atyp_hyp)<-"Ever had biopsy with atypical hyperplasia"

#Diagnosis variables
label(valdata$case)<-"BC case, anytime"
label(valdata$bca_invasive)<-"Invasive BC, anytime"
label(valdata$dxage_yr)<-"Age at diagnosis (y)"
label(valdata$case_5y)<-"BC case within 5y"
label(valdata$case_5yinv)<-"Invasive BC case within 5y"
label(valdata$case_10y)<-"BC case within 10y"
label(valdata$case_10yinv)<-"Invasive BC case within 10y"

table1<-table1(~ age+birth_year+race+height+bmi+bmi_cat4+agemenarche+parous+parity_new+ageflb+ocuse+meno_status2+hrt_use+pmh_type_any+
         alcohol_status+alcgm+smoking_status+smoking_amt+famhx_first+ever_biopsy+biopsy_num+bbd_history+atyp_hyp+
         case+bca_invasive+dxage_yr+case_5y+case_5yinv+case_10y+case_10yinv, data=valdata)
table1
#PLEASE INSERT FILEPATH AND THE APPROPRIATE ACRONYM 
write.table(table1, "FILEPATH/Table1_ACRONYM.csv", col.names=T, row.names=F, append=T, sep=',')

#Please send resulting table to kristen_brantley@dfci.harvard.edu