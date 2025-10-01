#TC datasetup for external cohorts
#Please pay careful attention to instructions at the end of this file when creating TC dataset
#############################################
# TYRER CUZICK
##############################################
rm(list=ls())

if(!require(tidyverse)){
  install.packages("tidyverse")
  suppress.warnings(library("dplyr"))
}
library(iCARE)
library(stringr)

#############
#Step 1. Create file in proper format.
###############
#output of risk file is:
#personal risk 5y
#population risk 5 y
#probability of brca1 gene
#probability of brca2 gene
#Model calibration datasets are already all clean

#Create new variables for tyrer formatting:
#Make sure all rounding for ages is floor
names(valdata)
drops<-c("V1","X")
tyrer_data<-valdata[,(! names(valdata) %in% drops),]
tyrer_data2<-valdata2[,(! names(valdata2) %in% drops),]

#Dataset #1: 
#Assumptions:
#Assume all with famhx=1 have a mother with BC unless other information was provided
#age at BC was 45 if <50 and 60 otherwise. If no information on fam1grBC50 variable, then assume age @ dx was 50

data_tyrer_assump1 <- tyrer_data%>%
  mutate(ty_id= gsub("_","",subject_id),
         ty_age=ifelse(!is.na(age), as.integer(floor(age)), -99),
         ty_agemenarche_cont=ifelse(!is.na(agemenarche), floor(agemenarche), -99),
         ty_parity=ifelse(!is.na(parity_new), parity_new, 2),
         ty_agefb=ifelse(!is.na(parity_new) & !is.na(ageflb), floor(ageflb), 
                         ifelse(is.na(parity_new)|parity_new==0, 0, -99)),
         ty_meno=ifelse(!is.na(meno_status) & meno_status==1, 2,
                        ifelse(!is.na(meno_status) & meno_status==2, 0,
                               ifelse(!is.na(meno_status) & meno_status==3, 1, 3))),
         ty_agemeno=ifelse(!is.na(meno_age) & ty_meno==2, floor(meno_age), -99),
         ty_heightm=ifelse(!is.na(height), round((height/100),digits=2), -99),
         ty_weightkg=ifelse(!is.na(weight), round(weight,digits=2), -99),
         ty_hyperplasia=ifelse(!is.na(bbd_history) & !is.na(bbd_type1) & bbd_history==1 & bbd_type1==2, 1, 0),
         ty_atyp_hyp=ifelse(! is.na(atyp_hyp) & atyp_hyp==1, 1, 0),
         ty_lcis=0,
         ty_blank=0,
         #individual had ovarian cancer- 0 for all b/c not collected
         ty_hist_ovarian=0,
         ty_age_ovarian=-99,
         ty_ashkenazi=ifelse(!is.na(ajancestry) & ajancestry==1, 1, 0),
         #hrt_use=c("NA"=NA, "Never use"=1, "Current use"=2, "Former use"=3, "Ever use (unk)"=4),
         #assume if former use or ever use they are previous user with >5 years since last use. 
         #assume never for those with missing
         ty_hrtuse1=ifelse( !is.na(hrt_use) & hrt_use==2, 3, 
                            ifelse(!is.na(hrt_use) & (hrt_use==3|hrt_use==4), 1, 0)),
         #assume if former use or ever use they are previous user with <5 years since last use. 
         #ty_hrtuse2=ifelse(hrt_use==2, 3, 
         #                 ifelse(hrt_use==3|hrt_use==4, 2, 0)),
         
         #pmh_type_any=c("NA"=NA, "Never"=1, "Past, unk type"=2, "Past, E only"=3, "Past, Combined"=4, "Current, unk type"=5, "Current, E only"=6, "Current, combined"=7),
         #NA, 1, 2, 4, 5, 7 should all be 1 for TC
         #Here based on above assumption all HRT use is >5 years in past- set hrttype to 1 for allpast users
         ty_hrttype=ifelse(!is.na(hrt_use) & ty_hrtuse1==3 & (pmh_type_any==3| pmh_type_any==6), 0, 1),
         ty_hrttype=ifelse(is.na(ty_hrttype),1,ty_hrttype),
         ty_durhrt=ifelse(!is.na(hrt_dur) & (ty_hrtuse1==3 | ty_hrtuse1==1 | ty_hrtuse1==2), round((hrt_dur/12), digits=1), 0),
         ty_intended=0,
         ty_timesincehrt=0,
         ty_genetictest=0,
         ty_genetictestfath=0,
         
         #If no data available set to 0
         #Note: If first degree relative w/ BC assume mother (need to alter this if they have info on brCancerMom, Sis, Dau, etc.)
         ty_mother_bc= ifelse(famhx_first==1 & !is.na(famhx_first) & (is.na(brcancermom)| (brcancermom==0)| (brcancermom==1)) &
                                (is.na(brcancerdad)| brcancerdad==0) &(is.na(brcancersis)| brcancersis==0) &
                                (is.na(brcancerdau)| brcancerdau==0), 1,0),
         #mother had bilateral bc? unk. for all
         ty_mother_bc_bilat=0,
         #mother ovarian cancer?
         ty_mother_ovarian=ifelse(ovcancermom==1 & !is.na(ovcancermom), 1, 0),
         #age mother had bc - Set to 45 if before age 50, 60 if after 50, otherwise unknown.
         ty_mother_age_bc=ifelse(ty_mother_bc==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                                 ifelse(ty_mother_bc==1 & fam1grbc50==0 & !is.na(fam1grbc50), 60,
                                        ifelse(ty_mother_bc==1 & is.na(fam1grbc50), 50, -99))),
         ty_mother_age_bilat_bc=-99,
         ty_mother_age_ovar=-99,
         ty_genetictestmoth=0,
         ty_number_sis=ifelse(!is.na(sisters) & sisters>0, sisters, 0),
         #once we collect updated information on famhx this will change
         #NOTE: need to figure out coding for # of sisters.
         #set sister age at BC dx as proband age 
         #only needed if sister >=1
         ty_sis_bc=ifelse(!is.na(brcancersis) & brcancersis>0, 1, 0),
         ty_sis_bilat=0,
         ty_sis_ovar=ifelse(!is.na(ovcancersis) & ovcancersis==1, 1, 0),
         ty_sis_age=ifelse(ty_sis_bc==1, age, -99),
         ty_sis_agebilat=-99,
         ty_sis_ageovar=-99,
         ty_sis_genetic=0,
         
         ty_patgran=0,
         ty_patgranbilat=0,
         ty_patgranovar=0,
         ty_patgranage=-99,
         ty_patgranagebilat=-99,
         ty_patgranageovar=-99,
         ty_patgrangenetic=0,
         ty_matgran=0,
         ty_matgranbilat=0,
         ty_matgranovar=0,
         ty_matgranage=-99,
         ty_matgranagebilat=-99,
         ty_matgranageovar=-99,
         ty_matgrangenetic=0,
         
         #for this will have to assume they have 0 paternal and maternal aunts since no info
         ty_pataunts=0,
         # ty_patauntbc=0,
         # ty_patauntbilat=0,
         # ty_patauntovar=0,
         # ty_patauntage=-99,
         # ty_patauntagebilat=-99,
         # ty_patauntageovar=-99,
         # ty_patauntgenetic=0,
         
         ty_mataunts=0,
         #ty_matauntbc=0,
         #ty_matauntbilat=0,
         #ty_matauntovar=0,
         #ty_matauntage=-99,
         #ty_matauntagebilat=-99,
         #ty_matauntageovar=-99,
         #ty_matauntgenetic=0,
         
         #for this will need to assume 0 daughters bc have no info
         ty_daughters=ifelse(daughters>0 & !is.na(daughters), daughters, 0),
         
         #only needed if daughter >=1
         ty_daughter_bc=ifelse(brcancerdau>0 & !is.na(brcancerdau), 1, 0),
         ty_daughterbilat=0,
         ty_daughterovar=0,
         #assume daughters age is proband age-30 y. 
         ty_daughterage=ifelse(ty_daughter_bc==1, age-30, -99),
         ty_daughteragebilat=-99,
         ty_daughterageovar=-99,
         ty_daughtergenetic=0,
         
         #Third stage famhx info
         #No information provided so input 7 zeros
         ty_noinfo1=0,
         ty_noinfo2=0,
         ty_noinfo3=0,
         ty_noinfo4=0,
         ty_noinfo5=0,
         ty_noinfo6=0,
         ty_noinfo7=0, 
         
         #Last fields- correspond with table 4
         
         ty_fatherbc=ifelse(brcancerdad==1 & !is.na(brcancerdad), 1, 0),
         ty_brotherbc=0,
         ty_density=0, #can update with cohorts with density info - assume 0 for now for comparisons
         ty_densityval=-99 , 
         ty_polyscore= 0)


###############
#ASSUMPTION #2
###############

#Assume all with famhx=1 is the maternal grandmother  unless other information was provided
#age at BC was 45 if <50 and 60 otherwise. If no information on fam1grBC50 variable, then assume age @ dx was 50

data_tyrer_assump2 <- tyrer_data %>%
  mutate(ty_id= sub("_","",subject_id),
         ty_age=ifelse(!is.na(age), floor(age), -99),
         ty_agemenarche_cont=ifelse(!is.na(agemenarche), floor(agemenarche), -99),
         ty_parity=ifelse(!is.na(parity_new), parity_new, 2),
         ty_agefb=ifelse(!is.na(parity_new) & !is.na(ageflb), floor(ageflb), 
                         ifelse(is.na(parity_new)|parity_new==0, 0, -99)),
         ty_meno=ifelse(!is.na(meno_status) & meno_status==1, 2,
                        ifelse(!is.na(meno_status) & meno_status==2, 0,
                               ifelse(!is.na(meno_status) & meno_status==3, 1, 3))),
         ty_agemeno=ifelse(!is.na(meno_age) & ty_meno==2, floor(meno_age), -99),
         ty_heightm=ifelse(!is.na(height), round((height/100),digits=2), -99),
         ty_weightkg=ifelse(!is.na(weight), round(weight,digits=2), -99),
         ty_hyperplasia=ifelse(!is.na(bbd_history) & !is.na(bbd_type1) & bbd_history==1 & bbd_type1==2, 1, 0),
         ty_atyp_hyp=ifelse(! is.na(atyp_hyp) & atyp_hyp==1, 1, 0),
         ty_lcis=0,
         ty_blank=0,
         #individual had ovarian cancer- 0 for all b/c not collected
         ty_hist_ovarian=0,
         ty_age_ovarian=-99,
         ty_ashkenazi=ifelse(!is.na(ajancestry) & ajancestry==1, 1, 0),
         #hrt_use=c("NA"=NA, "Never use"=1, "Current use"=2, "Former use"=3, "Ever use (unk)"=4),
         #assume if former use or ever use they are previous user with >5 years since last use. 
         #assume never for those with missing
         ty_hrtuse1=ifelse( !is.na(hrt_use) & hrt_use==2, 3, 
                            ifelse(!is.na(hrt_use) & (hrt_use==3|hrt_use==4), 1, 0)),
         #assume if former use or ever use they are previous user with <5 years since last use. 
         #ty_hrtuse2=ifelse(hrt_use==2, 3, 
         #                 ifelse(hrt_use==3|hrt_use==4, 2, 0)),
         
         #pmh_type_any=c("NA"=NA, "Never"=1, "Past, unk type"=2, "Past, E only"=3, "Past, Combined"=4, "Current, unk type"=5, "Current, E only"=6, "Current, combined"=7),
         #NA, 1, 2, 4, 5, 7 should all be 1 for TC
         #Here based on above assumption all HRT use is >5 years in past- set hrttype to 1 for allpast users
         ty_hrttype=ifelse(!is.na(hrt_use) & ty_hrtuse1==3 & (pmh_type_any==3| pmh_type_any==6), 0, 1),
         ty_hrttype=ifelse(is.na(ty_hrttype),1,ty_hrttype),
         ty_durhrt=ifelse(!is.na(hrt_dur) & (ty_hrtuse1==3 | ty_hrtuse1==1 | ty_hrtuse1==2), round((hrt_dur/12), digits=1), 0),
         ty_durhrt=ifelse(ty_durhrt<1, round(ty_durhrt, digits=0), ty_durhrt),
         ty_intended=0,
         ty_timesincehrt=0,
         ty_genetictest=0,
         ty_genetictestfath=0,
         
         #If no data available set to 0
         #mother with bc?
         #Assumption #2: If first degree relative but no info, here we set as grandmother as the family member
         ty_mother_bc= ifelse(brcancermom==1 & !is.na(brcancermom), 1, 0),
         #mother had bilateral bc? unk. for all
         ty_mother_bc_bilat=0,
         #mother ovarian cancer?
         ty_mother_ovarian=ifelse(ovcancermom==1 & !is.na(ovcancermom), 1, 0),
         #age mother had bc - unknown for all
         ty_mother_age_bc=ifelse(ty_mother_bc==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                                 ifelse(ty_mother_bc==1 & fam1grbc50==0 & !is.na(fam1grbc50),60,
                                        ifelse(ty_mother_bc==1 & is.na(fam1grbc50), 50, -99))),
         ty_mother_age_bilat_bc=-99,
         ty_mother_age_ovar=-99,
         ty_genetictestmoth=0,
         ty_number_sis=ifelse(!is.na(sisters) & sisters>0, sisters, 0),
         #once we collect updated information on famhx this will change
         #NOTE: need to figure out coding for # of sisters.
         #only needed if sister >=1
         #Need to set all sister study to 1 for sister,not assume mother BC if famhx=1
         ty_sis_bc=ifelse(!is.na(brcancersis) & brcancersis>0, 1, 0),
         ty_sis_bilat=0,
         ty_sis_ovar=ifelse(!is.na(ovcancersis) & ovcancersis==1, 1, 0),
         #set to age of proband for sister
         ty_sis_age=ifelse(ty_sis_bc==1, age, -99),
         ty_sis_agebilat=-99,
         ty_sis_ageovar=-99,
         ty_sis_genetic=0,
         
         ty_patgran=0,
         ty_patgranbilat=0,
         ty_patgranovar=0,
         ty_patgranage=-99,
         ty_patgranagebilat=-99,
         ty_patgranageovar=-99,
         ty_patgrangenetic=0,
         ty_matgran=ifelse(famhx_first==1 & !is.na(famhx_first) & (is.na(brcancermom)| (brcancermom==0)) &
                             (is.na(brcancerdad)| brcancerdad==0) &(is.na(brcancersis)| brcancersis==0) &
                             (is.na(brcancerdau)| brcancerdau==0), 1,0),
         ty_matgranbilat=0,
         ty_matgranovar=0,
         ty_matgranage=ifelse(ty_matgran==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                              ifelse(ty_matgran==1 & fam1grbc50==0 & !is.na(fam1grbc50), 60,
                                     ifelse(ty_matgran==1 & is.na(fam1grbc50), 50, -99))),
         ty_matgranagebilat=-99,
         ty_matgranageovar=-99,
         ty_matgrangenetic=0,
         
         #for this will have to assume they have 0 paternal and maternal aunts since no info
         ty_pataunts=0,
         # ty_patauntbc=0,
         # ty_patauntbilat=0,
         # ty_patauntovar=0,
         # ty_patauntage=-99,
         # ty_patauntagebilat=-99,
         # ty_patauntageovar=-99,
         # ty_patauntgenetic=0,
         
         ty_mataunts=0,
         #ty_matauntbc=0,
         #ty_matauntbilat=0,
         #ty_matauntovar=0,
         #ty_matauntage=-99,
         #ty_matauntagebilat=-99,
         #ty_matauntageovar=-99,
         #ty_matauntgenetic=0,
         
         #for this will need to assume 0 daughters bc have no info
         ty_daughters=ifelse(daughters>0 & !is.na(daughters), daughters, 0),
         
         #only needed if daughter >=1
         ty_daughter_bc=ifelse(brcancerdau>0 & !is.na(brcancerdau), 1, 0),
         ty_daughterbilat=0,
         ty_daughterovar=0,
         ty_daughterage=ifelse(ty_daughter_bc==1, age-30, -99),
         ty_daughteragebilat=-99,
         ty_daughterageovar=-99,
         ty_daughtergenetic=0,
         
         #Third stage famhx info
         #No information provided so input 7 zeros
         ty_noinfo1=0,
         ty_noinfo2=0,
         ty_noinfo3=0,
         ty_noinfo4=0,
         ty_noinfo5=0,
         ty_noinfo6=0,
         ty_noinfo7=0, 
         
         #Last fields- correspond with table 4
         
         ty_fatherbc=ifelse(brcancerdad==1 & !is.na(brcancerdad), 1, 0),
         ty_brotherbc=0,
         ty_density=0, #can update with cohorts with density info - assume 0 for now for comparisons
         ty_densityval=-99 , 
         ty_polyscore= 0)

#2 iterations: assume mom, assume maternal grandma

names<-colnames(data_tyrer_assump1)
names_tyrer<-names[grepl("^ty_", names)]
names_tyrer

#For those assuming mother is with BC (all except SIS, BWHS, MMHS) will take out variable names for sisters or daughters
names_tyrer_ex_sis<-c("ty_sis_bc", "ty_sis_bilat","ty_sis_ovar","ty_sis_age","ty_sis_agebilat",
                      "ty_sis_ageovar","ty_sis_genetic")

names_tyrer_ex_dau<-c("ty_daughter_bc", "ty_daughterbilat", "ty_daughterovar", "ty_daughterage", 
                      "ty_daughteragebilat","ty_daughterageovar","ty_daughtergenetic")

names_tyrer2<-names_tyrer[! names_tyrer %in% names_tyrer_ex_sis]
names_tyrer2
names_tyrer3<-names_tyrer[! names_tyrer %in% names_tyrer_ex_dau]
names_tyrer3
names_tyrer_ex_sisdau<-c(names_tyrer_ex_sis, names_tyrer_ex_dau)
names_tyrer4<-names_tyrer[! names_tyrer %in% names_tyrer_ex_sisdau]

# #Need to make sure the uploaded list is specific to each id based on whether they have daughters or sisters. 


tyrer1_modeldata_assump1<-data_tyrer_assump1[,(names(data_tyrer_assump1) %in% names_tyrer4)]
tyrer1_modeldata_assump2<-data_tyrer_assump2[,(names(data_tyrer_assump2) %in% names_tyrer4)]

#########################################################################
###INSTRUCTIONS: PLEASE READ.
## YOU NEED TO SEE WHAT SIZE YOUR DATASET IS BEFORE THIS NEXT PART
## The calculator only seems to handle up to 25,000 observations at a time, so you will need to break down, as demonstrated below.
## You will need to change the #s provided in []s . This is an example of a dataset with 101256 individuals.
#This creates datasets for 1-5y
#You also need to provide FILEPATHs for all of these - within the same folder everything else is stored in.
#################################################

#Find # of observations
dim(tyrer1_modeldata_assump1)

#Seems to handle up to 25,000 at a time - more than that and it runs out of memory
tyrer1_assump1_1<-tyrer1_modeldata_assump1[1:25000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump1_1.txt")
write.table(tyrer1_assump1_1,"FILEPATH.data_1to5y_assump1_1.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer1_assump1_2<-tyrer1_modeldata_assump1[25001:50000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump1_2.txt")
write.table(tyrer1_assump1_2,"FILEPATH/data_1to5y_assump1_2.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer1_assump1_3<-tyrer1_modeldata_assump1[50001:75000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump1_3.txt")
write.table(tyrer1_assump1_3,"FILEPATH/data_1to5y_assump1_3.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer1_assump1_4<-tyrer1_modeldata_assump1[75001:100000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump1_4.txt")
write.table(tyrer1_assump1_3,"FILEPATH/data_1to5y_assump1_4.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

##READ THIS! Make sure you change the number given after "v8" here so that it matches the # of observations in your dataset

tyrer1_assump1_5<-tyrer1_modeldata_assump1[100001:101256,]
writeLines(c("v8","1256"), "FILEPATH/data_1to5y_assump1_5.txt")
write.table(tyrer1_assump1_3,"FILEPATH/data_1to5y_assump1_5.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 


######
#Assumption 2
tyrer1_assump2_1<-tyrer1_modeldata_assump2[1:25000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump2_1.txt")
write.table(tyrer1_assump2_1,"FILEPATH.data_1to5y_assump2_1.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer1_assump2_2<-tyrer1_modeldata_assump2[25001:50000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump2_2.txt")
write.table(tyrer1_assump2_2,"FILEPATH/data_1to5y_assump2_2.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer1_assump2_3<-tyrer1_modeldata_assump2[50001:75000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump2_3.txt")
write.table(tyrer1_assump2_3,"FILEPATH/data_1to5y_assump2_3.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer1_assump2_4<-tyrer1_modeldata_assump2[75001:100000,]
writeLines(c("v8","25000"), "FILEPATH/data_1to5y_assump2_4.txt")
write.table(tyrer1_assump2_3,"FILEPATH/data_1to5y_assump2_4.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

##READ THIS! Make sure you change the number given after "v8" here so that it matches the # of observations in your dataset

tyrer1_assump2_5<-tyrer1_modeldata_assump2[100001:101256,]
writeLines(c("v8","1256"), "FILEPATH/data_1to5y_assump2_5.txt")
write.table(tyrer1_assump2_3,"FILEPATH/data_1to5y_assump2_5.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 


#############################
# 2-6y
###############################
#Dataset #1: 
#Assumptions:
#Assume all with famhx=1 have a mother with BC unless other information was provided
#age at BC was 45 if <50 and 60 otherwise. If no information on fam1grBC50 variable, then assume age @ dx was 50

data_tyrer_assump1 <- tyrer_data2 %>%
  mutate(ty_id= gsub("_","",subject_id),
         ty_age=ifelse(!is.na(age), as.integer(floor(age)), -99),
         ty_agemenarche_cont=ifelse(!is.na(agemenarche), floor(agemenarche), -99),
         ty_parity=ifelse(!is.na(parity_new), parity_new, 2),
         ty_agefb=ifelse(!is.na(parity_new) & !is.na(ageflb), floor(ageflb), 
                         ifelse(is.na(parity_new)|parity_new==0, 0, -99)),
         ty_meno=ifelse(!is.na(meno_status) & meno_status==1, 2,
                        ifelse(!is.na(meno_status) & meno_status==2, 0,
                               ifelse(!is.na(meno_status) & meno_status==3, 1, 3))),
         ty_agemeno=ifelse(!is.na(meno_age) & ty_meno==2, floor(meno_age), -99),
         ty_heightm=ifelse(!is.na(height), round((height/100),digits=2), -99),
         ty_weightkg=ifelse(!is.na(weight), round(weight,digits=2), -99),
         ty_hyperplasia=ifelse(!is.na(bbd_history) & !is.na(bbd_type1) & bbd_history==1 & bbd_type1==2, 1, 0),
         ty_atyp_hyp=ifelse(! is.na(atyp_hyp) & atyp_hyp==1, 1, 0),
         ty_lcis=0,
         ty_blank=0,
         #individual had ovarian cancer- 0 for all b/c not collected
         ty_hist_ovarian=0,
         ty_age_ovarian=-99,
         ty_ashkenazi=ifelse(!is.na(ajancestry) & ajancestry==1, 1, 0),
         #hrt_use=c("NA"=NA, "Never use"=1, "Current use"=2, "Former use"=3, "Ever use (unk)"=4),
         #assume if former use or ever use they are previous user with >5 years since last use. 
         #assume never for those with missing
         ty_hrtuse1=ifelse( !is.na(hrt_use) & hrt_use==2, 3, 
                            ifelse(!is.na(hrt_use) & (hrt_use==3|hrt_use==4), 1, 0)),
         #assume if former use or ever use they are previous user with <5 years since last use. 
         #ty_hrtuse2=ifelse(hrt_use==2, 3, 
         #                 ifelse(hrt_use==3|hrt_use==4, 2, 0)),
         
         #pmh_type_any=c("NA"=NA, "Never"=1, "Past, unk type"=2, "Past, E only"=3, "Past, Combined"=4, "Current, unk type"=5, "Current, E only"=6, "Current, combined"=7),
         #NA, 1, 2, 4, 5, 7 should all be 1 for TC
         #Here based on above assumption all HRT use is >5 years in past- set hrttype to 1 for allpast users
         ty_hrttype=ifelse(!is.na(hrt_use) & ty_hrtuse1==3 & (pmh_type_any==3| pmh_type_any==6), 0, 1),
         ty_hrttype=ifelse(is.na(ty_hrttype),1,ty_hrttype),
         ty_durhrt=ifelse(!is.na(hrt_dur) & (ty_hrtuse1==3 | ty_hrtuse1==1 | ty_hrtuse1==2), round((hrt_dur/12), digits=1), 0),
         ty_intended=0,
         ty_timesincehrt=0,
         ty_genetictest=0,
         ty_genetictestfath=0,
         
         #If no data available set to 0
         #Note: If first degree relative w/ BC assume mother (need to alter this if they have info on brCancerMom, Sis, Dau, etc.)
         ty_mother_bc= ifelse(famhx_first==1 & !is.na(famhx_first) & (is.na(brcancermom)| (brcancermom==0)| (brcancermom==1)) &
                                (is.na(brcancerdad)| brcancerdad==0) &(is.na(brcancersis)| brcancersis==0) &
                                (is.na(brcancerdau)| brcancerdau==0), 1,0),
         #mother had bilateral bc? unk. for all
         ty_mother_bc_bilat=0,
         #mother ovarian cancer?
         ty_mother_ovarian=ifelse(ovcancermom==1 & !is.na(ovcancermom), 1, 0),
         #age mother had bc - Set to 45 if before age 50, 60 if after 50, otherwise unknown.
         ty_mother_age_bc=ifelse(ty_mother_bc==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                                 ifelse(ty_mother_bc==1 & fam1grbc50==0 & !is.na(fam1grbc50), 60,
                                        ifelse(ty_mother_bc==1 & is.na(fam1grbc50), 50, -99))),
         ty_mother_age_bilat_bc=-99,
         ty_mother_age_ovar=-99,
         ty_genetictestmoth=0,
         ty_number_sis=ifelse(!is.na(sisters) & sisters>0, sisters, 0),
         #once we collect updated information on famhx this will change
         #NOTE: need to figure out coding for # of sisters.
         #set sister age at BC dx as proband age 
         #only needed if sister >=1
         ty_sis_bc=ifelse(!is.na(brcancersis) & brcancersis>0, 1, 0),
         ty_sis_bilat=0,
         ty_sis_ovar=ifelse(!is.na(ovcancersis) & ovcancersis==1, 1, 0),
         ty_sis_age=ifelse(ty_sis_bc==1, age, -99),
         ty_sis_agebilat=-99,
         ty_sis_ageovar=-99,
         ty_sis_genetic=0,
         
         ty_patgran=0,
         ty_patgranbilat=0,
         ty_patgranovar=0,
         ty_patgranage=-99,
         ty_patgranagebilat=-99,
         ty_patgranageovar=-99,
         ty_patgrangenetic=0,
         ty_matgran=0,
         ty_matgranbilat=0,
         ty_matgranovar=0,
         ty_matgranage=-99,
         ty_matgranagebilat=-99,
         ty_matgranageovar=-99,
         ty_matgrangenetic=0,
         
         #for this will have to assume they have 0 paternal and maternal aunts since no info
         ty_pataunts=0,
         # ty_patauntbc=0,
         # ty_patauntbilat=0,
         # ty_patauntovar=0,
         # ty_patauntage=-99,
         # ty_patauntagebilat=-99,
         # ty_patauntageovar=-99,
         # ty_patauntgenetic=0,
         
         ty_mataunts=0,
         #ty_matauntbc=0,
         #ty_matauntbilat=0,
         #ty_matauntovar=0,
         #ty_matauntage=-99,
         #ty_matauntagebilat=-99,
         #ty_matauntageovar=-99,
         #ty_matauntgenetic=0,
         
         #for this will need to assume 0 daughters bc have no info
         ty_daughters=ifelse(daughters>0 & !is.na(daughters), daughters, 0),
         
         #only needed if daughter >=1
         ty_daughter_bc=ifelse(brcancerdau>0 & !is.na(brcancerdau), 1, 0),
         ty_daughterbilat=0,
         ty_daughterovar=0,
         #assume daughters age is proband age-30 y. 
         ty_daughterage=ifelse(ty_daughter_bc==1, age-30, -99),
         ty_daughteragebilat=-99,
         ty_daughterageovar=-99,
         ty_daughtergenetic=0,
         
         #Third stage famhx info
         #No information provided so input 7 zeros
         ty_noinfo1=0,
         ty_noinfo2=0,
         ty_noinfo3=0,
         ty_noinfo4=0,
         ty_noinfo5=0,
         ty_noinfo6=0,
         ty_noinfo7=0, 
         
         #Last fields- correspond with table 4
         
         ty_fatherbc=ifelse(brcancerdad==1 & !is.na(brcancerdad), 1, 0),
         ty_brotherbc=0,
         ty_density=0, #can update with cohorts with density info - assume 0 for now for comparisons
         ty_densityval=-99 , 
         ty_polyscore= 0)


###############
#ASSUMPTION #2
###############

#Assume all with famhx=1 is the maternal grandmother  unless other information was provided
#age at BC was 45 if <50 and 60 otherwise. If no information on fam1grBC50 variable, then assume age @ dx was 50

data_tyrer_assump2 <- tyrer_data2 %>%
  mutate(ty_id= sub("_","",subject_id),
         ty_age=ifelse(!is.na(age), floor(age), -99),
         ty_agemenarche_cont=ifelse(!is.na(agemenarche), floor(agemenarche), -99),
         ty_parity=ifelse(!is.na(parity_new), parity_new, 2),
         ty_agefb=ifelse(!is.na(parity_new) & !is.na(ageflb), floor(ageflb), 
                         ifelse(is.na(parity_new)|parity_new==0, 0, -99)),
         ty_meno=ifelse(!is.na(meno_status) & meno_status==1, 2,
                        ifelse(!is.na(meno_status) & meno_status==2, 0,
                               ifelse(!is.na(meno_status) & meno_status==3, 1, 3))),
         ty_agemeno=ifelse(!is.na(meno_age) & ty_meno==2, floor(meno_age), -99),
         ty_heightm=ifelse(!is.na(height), round((height/100),digits=2), -99),
         ty_weightkg=ifelse(!is.na(weight), round(weight,digits=2), -99),
         ty_hyperplasia=ifelse(!is.na(bbd_history) & !is.na(bbd_type1) & bbd_history==1 & bbd_type1==2, 1, 0),
         ty_atyp_hyp=ifelse(! is.na(atyp_hyp) & atyp_hyp==1, 1, 0),
         ty_lcis=0,
         ty_blank=0,
         #individual had ovarian cancer- 0 for all b/c not collected
         ty_hist_ovarian=0,
         ty_age_ovarian=-99,
         ty_ashkenazi=ifelse(!is.na(ajancestry) & ajancestry==1, 1, 0),
         #hrt_use=c("NA"=NA, "Never use"=1, "Current use"=2, "Former use"=3, "Ever use (unk)"=4),
         #assume if former use or ever use they are previous user with >5 years since last use. 
         #assume never for those with missing
         ty_hrtuse1=ifelse( !is.na(hrt_use) & hrt_use==2, 3, 
                            ifelse(!is.na(hrt_use) & (hrt_use==3|hrt_use==4), 1, 0)),
         #assume if former use or ever use they are previous user with <5 years since last use. 
         #ty_hrtuse2=ifelse(hrt_use==2, 3, 
         #                 ifelse(hrt_use==3|hrt_use==4, 2, 0)),
         
         #pmh_type_any=c("NA"=NA, "Never"=1, "Past, unk type"=2, "Past, E only"=3, "Past, Combined"=4, "Current, unk type"=5, "Current, E only"=6, "Current, combined"=7),
         #NA, 1, 2, 4, 5, 7 should all be 1 for TC
         #Here based on above assumption all HRT use is >5 years in past- set hrttype to 1 for allpast users
         ty_hrttype=ifelse(!is.na(hrt_use) & ty_hrtuse1==3 & (pmh_type_any==3| pmh_type_any==6), 0, 1),
         ty_hrttype=ifelse(is.na(ty_hrttype),1,ty_hrttype),
         ty_durhrt=ifelse(!is.na(hrt_dur) & (ty_hrtuse1==3 | ty_hrtuse1==1 | ty_hrtuse1==2), round((hrt_dur/12), digits=1), 0),
         ty_durhrt=ifelse(ty_durhrt<1, round(ty_durhrt, digits=0), ty_durhrt),
         ty_intended=0,
         ty_timesincehrt=0,
         ty_genetictest=0,
         ty_genetictestfath=0,
         
         #If no data available set to 0
         #mother with bc?
         #Assumption #2: If first degree relative but no info, here we set as grandmother as the family member
         ty_mother_bc= ifelse(brcancermom==1 & !is.na(brcancermom), 1, 0),
         #mother had bilateral bc? unk. for all
         ty_mother_bc_bilat=0,
         #mother ovarian cancer?
         ty_mother_ovarian=ifelse(ovcancermom==1 & !is.na(ovcancermom), 1, 0),
         #age mother had bc - unknown for all
         ty_mother_age_bc=ifelse(ty_mother_bc==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                                 ifelse(ty_mother_bc==1 & fam1grbc50==0 & !is.na(fam1grbc50),60,
                                        ifelse(ty_mother_bc==1 & is.na(fam1grbc50), 50, -99))),
         ty_mother_age_bilat_bc=-99,
         ty_mother_age_ovar=-99,
         ty_genetictestmoth=0,
         ty_number_sis=ifelse(!is.na(sisters) & sisters>0, sisters, 0),
         #once we collect updated information on famhx this will change
         #NOTE: need to figure out coding for # of sisters.
         #only needed if sister >=1
         #Need to set all sister study to 1 for sister,not assume mother BC if famhx=1
         ty_sis_bc=ifelse(!is.na(brcancersis) & brcancersis>0, 1, 0),
         ty_sis_bilat=0,
         ty_sis_ovar=ifelse(!is.na(ovcancersis) & ovcancersis==1, 1, 0),
         #set to age of proband for sister
         ty_sis_age=ifelse(ty_sis_bc==1, age, -99),
         ty_sis_agebilat=-99,
         ty_sis_ageovar=-99,
         ty_sis_genetic=0,
         
         ty_patgran=0,
         ty_patgranbilat=0,
         ty_patgranovar=0,
         ty_patgranage=-99,
         ty_patgranagebilat=-99,
         ty_patgranageovar=-99,
         ty_patgrangenetic=0,
         ty_matgran=ifelse(famhx_first==1 & !is.na(famhx_first) & (is.na(brcancermom)| (brcancermom==0)) &
                             (is.na(brcancerdad)| brcancerdad==0) &(is.na(brcancersis)| brcancersis==0) &
                             (is.na(brcancerdau)| brcancerdau==0), 1,0),
         ty_matgranbilat=0,
         ty_matgranovar=0,
         ty_matgranage=ifelse(ty_matgran==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                              ifelse(ty_matgran==1 & fam1grbc50==0 & !is.na(fam1grbc50), 60,
                                     ifelse(ty_matgran==1 & is.na(fam1grbc50), 50, -99))),
         ty_matgranagebilat=-99,
         ty_matgranageovar=-99,
         ty_matgrangenetic=0,
         
         #for this will have to assume they have 0 paternal and maternal aunts since no info
         ty_pataunts=0,
         # ty_patauntbc=0,
         # ty_patauntbilat=0,
         # ty_patauntovar=0,
         # ty_patauntage=-99,
         # ty_patauntagebilat=-99,
         # ty_patauntageovar=-99,
         # ty_patauntgenetic=0,
         
         ty_mataunts=0,
         #ty_matauntbc=0,
         #ty_matauntbilat=0,
         #ty_matauntovar=0,
         #ty_matauntage=-99,
         #ty_matauntagebilat=-99,
         #ty_matauntageovar=-99,
         #ty_matauntgenetic=0,
         
         #for this will need to assume 0 daughters bc have no info
         ty_daughters=ifelse(daughters>0 & !is.na(daughters), daughters, 0),
         
         #only needed if daughter >=1
         ty_daughter_bc=ifelse(brcancerdau>0 & !is.na(brcancerdau), 1, 0),
         ty_daughterbilat=0,
         ty_daughterovar=0,
         ty_daughterage=ifelse(ty_daughter_bc==1, age-30, -99),
         ty_daughteragebilat=-99,
         ty_daughterageovar=-99,
         ty_daughtergenetic=0,
         
         #Third stage famhx info
         #No information provided so input 7 zeros
         ty_noinfo1=0,
         ty_noinfo2=0,
         ty_noinfo3=0,
         ty_noinfo4=0,
         ty_noinfo5=0,
         ty_noinfo6=0,
         ty_noinfo7=0, 
         
         #Last fields- correspond with table 4
         
         ty_fatherbc=ifelse(brcancerdad==1 & !is.na(brcancerdad), 1, 0),
         ty_brotherbc=0,
         ty_density=0, #can update with cohorts with density info - assume 0 for now for comparisons
         ty_densityval=-99 , 
         ty_polyscore= 0)

#2 iterations: assume mom, assume maternal grandma

names<-colnames(data_tyrer_assump1)
names_tyrer<-names[grepl("^ty_", names)]
names_tyrer

#For those assuming mother is with BC (all except SIS, BWHS, MMHS) will take out variable names for sisters or daughters
names_tyrer_ex_sis<-c("ty_sis_bc", "ty_sis_bilat","ty_sis_ovar","ty_sis_age","ty_sis_agebilat",
                      "ty_sis_ageovar","ty_sis_genetic")

names_tyrer_ex_dau<-c("ty_daughter_bc", "ty_daughterbilat", "ty_daughterovar", "ty_daughterage", 
                      "ty_daughteragebilat","ty_daughterageovar","ty_daughtergenetic")

names_tyrer2<-names_tyrer[! names_tyrer %in% names_tyrer_ex_sis]
names_tyrer2
names_tyrer3<-names_tyrer[! names_tyrer %in% names_tyrer_ex_dau]
names_tyrer3
names_tyrer_ex_sisdau<-c(names_tyrer_ex_sis, names_tyrer_ex_dau)
names_tyrer4<-names_tyrer[! names_tyrer %in% names_tyrer_ex_sisdau]

# #Need to make sure the uploaded list is specific to each id based on whether they have daughters or sisters. 

tyrer2_modeldata_assump1<-data_tyrer_assump1[,(names(data_tyrer_assump1) %in% names_tyrer4)]
tyrer2_modeldata_assump2<-data_tyrer_assump2[,(names(data_tyrer_assump2) %in% names_tyrer4)]

#The size of your dataset will have changed. Please re-run dim to find the new # for the 2-6 y follow up.
dim(tyrer2_modeldata_assump1)

tyrer2_assump1_1<-tyrer2_modeldata_assump1[1:25000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump1_1.txt")
write.table(tyrer2_assump1_1,"FILEPATH.data_2to6y_assump1_1.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer2_assump1_2<-tyrer2_modeldata_assump1[25001:50000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump1_2.txt")
write.table(tyrer2_assump1_2,"FILEPATH/data_2to6y_assump1_2.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer2_assump1_3<-tyrer2_modeldata_assump1[50001:75000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump1_3.txt")
write.table(tyrer2_assump1_3,"FILEPATH/data_2to6y_assump1_3.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer2_assump1_4<-tyrer2_modeldata_assump1[75001:100000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump1_4.txt")
write.table(tyrer2_assump1_3,"FILEPATH/data_2to6y_assump1_4.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

##READ THIS! Make sure you change the number given after "v8" here so that it matches the # of observations in your dataset

tyrer2_assump1_5<-tyrer2_modeldata_assump1[100001:101256,]
writeLines(c("v8","1256"), "FILEPATH/data_2to6y_assump1_5.txt")
write.table(tyrer2_assump1_3,"FILEPATH/data_2to6y_assump1_5.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 


######
#Assumption 2
tyrer2_assump2_1<-tyrer2_modeldata_assump2[1:25000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump2_1.txt")
write.table(tyrer2_assump2_1,"FILEPATH.data_2to6y_assump2_1.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer2_assump2_2<-tyrer2_modeldata_assump2[25001:50000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump2_2.txt")
write.table(tyrer2_assump2_2,"FILEPATH/data_2to6y_assump2_2.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer2_assump2_3<-tyrer2_modeldata_assump2[50001:75000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump2_3.txt")
write.table(tyrer2_assump2_3,"FILEPATH/data_2to6y_assump2_3.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

tyrer2_assump2_4<-tyrer2_modeldata_assump2[75001:100000,]
writeLines(c("v8","25000"), "FILEPATH/data_2to6y_assump2_4.txt")
write.table(tyrer2_assump2_3,"FILEPATH/data_2to6y_assump2_4.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 

##READ THIS! Make sure you change the number given after "v8" here so that it matches the # of observations in your dataset

tyrer2_assump2_5<-tyrer2_modeldata_assump2[100001:101256,]
writeLines(c("v8","1256"), "FILEPATH/data_2to6y_assump2_5.txt")
write.table(tyrer2_assump2_3,"FILEPATH/data_2to6y_assump2_5.txt", col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") 
