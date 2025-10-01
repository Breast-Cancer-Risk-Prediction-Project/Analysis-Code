
#############################################
# TYRER CUZICK
##############################################
rm(list=ls())
########################################################
#Step 1. Read in box directory
########################################################
#Changed to invasive only - 5.22.24
install.packages("tidyverse")
install.packages("palmerpenguins")
install.packages(c("boxr", "base", "usethis"))
library(dplyr)
library(boxr)
library(purrr)
library(tidyr)

#If you haven't already go ahead and source this code:
source("C:/Users/krist/Dropbox/BCRPProject/Code/Model Calibration/Internal Code/Create2datasets.R")

#Now you should have a dataset called "valdata".

#Step 1. Create file in proper format.
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

library(stringr)


#Dataset 1-5 y (Repeat all below steps for dataset 2 (2-6y, tyrer_data2))

#Assumption 1:
#Assume all with famhx=1 have a mother with BC unless other information was provided
#age at BC was 45 if <50 and 60 otherwise. If no information on fam1grBC50 variable, then assume age @ dx was 50


#No. of sisters and daughters can only go up to 5 so set sisters_new to no. of sisters and if >5 set to 5
tyrer_data<-tyrer_data%>%
  rowwise() %>%
  #var for missing brcancer status when sisters are reported
  mutate(sister_nabc=ifelse(sisters>0 & is.na(brcancersis), sisters, 0),
         brcancersis_nomiss=ifelse(!is.na(brcancersis), brcancersis, 0),
         sister_nov=ifelse(sisters>0 & is.na(ovcancersis), sisters, 0),
         ovcancersis_nomiss=ifelse(!is.na(ovcancersis), ovcancersis, 0),
         
         daughter_nabc=ifelse(daughters>0 & is.na(brcancerdau), daughters,0),
         brcancerdau_nomiss=ifelse(!is.na(brcancerdau), brcancerdau,0),
         daughter_nov=ifelse(daughters>0 & is.na(daughters), daughters,0),
         ovcancerdau_nomiss=ifelse(!is.na(ovcancerdau), daughters, 0),
         sisters_new=ifelse(is.na(sisters), 0 , 
                            ifelse(!is.na(sisters) & sisters<=5, sisters,
                                   ifelse(sisters>5, 5))),
         daughters_new=ifelse(is.na(daughters),0,
                              ifelse(!is.na(daughters) & daughters<=5, daughters,
                                     ifelse(daughters>5, 5))))


#head(tyrer_data[1:20, c("sisters","brcancersis","brcancersis_nomiss","sister_nabc")])

#dummy columns - make unique to each individual - max # of sisters recorded
#Will have to assume that if sisters N is missing, then sisters=0 - use sisters_new

tyrer_data <- tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sister_bc = list({
      n <- sisters_new
      n_bc <-  brcancersis_nomiss
      if (n == 0) {
        rep(NA, sisters_new+1)  #if n=0 for sisters, then no information on sister w/ bc etc., so set to NA -will delete column later
      } else {
        #for sister_bc add a 1 for each sister that was=1
        #add a 0 for each sister that exists that did not have BC (n-n_bc)
        c(rep(1, n_bc), rep(0, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_sister_bc, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_bc", seq_along(.), ""), starts_with("ty_sister_bc"))

#check - this works!! 
#tyrer_data[which (tyrer_data$sisters_new==2),][c(1:30),c("sisters","brcancersis","ty_sister_bc1","ty_sister_bc2")]
#Repeat for other TC variables related to number of sisters

tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sister_bilat = list({
      n <- sisters_new
      n_bc <-  brcancersis_nomiss
      if (n == 0) {
        rep(NA, sisters_new+1) #Add +1 so variable NA gets created if 0 sisters - will delete
      } else {
        #for sister_bc add a 1 for each sister that was=1
        #add a 0 for each sister that exists that did not have BC (n-n_bc)
        #add NA for NA info on sisters with bc
        #add NA for all the sisters that are not included per individual (i.e., if they have 6 sisters, 6 dummy vars get NA)
        c(rep(0, n_bc), rep(0, n - n_bc))
        
      }
    })
  ) %>%
  unnest_wider(ty_sister_bilat, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_bilat", seq_along(.), ""), starts_with("ty_sister_bilat"))


tyrer_data<-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sister_age = list({
      n <- sisters_new
      n_bc <-  brcancersis_nomiss
      if (n == 0) {
        rep(NA, sisters_new+1) #Add +1 so variable gets created if 0 sisters
      } 
      else {
        
        #in BCRPP, sisters age not asked- leave as -99 for assumption 1 (assuming if famhx first=1 then it was mother)
        c(rep(-99, n_bc), rep(-99, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_sister_age, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_age", seq_along(.), ""), starts_with("ty_sister_age"))


#variable not collected in bcrpp
tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sister_bilatage = list({
      n <- sisters_new
      n_bc <-  brcancersis_nomiss
      if (n == 0 ) {
        rep(NA, sisters_new+1) #Add +1 so variable gets created if 0 sisters
      } 
      else {
        #for sister_bc add a 1 for each sister that was=1
        #add a 0 for each sister that exists that did not have BC (n-n_bc)
        #add NA for NA info on sisters with bc
        #add NA for all the sisters that are not included per individual (i.e., if they have 6 sisters, 6 dummy vars get NA)
        c(rep(-99, n_bc), rep(-99, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_sister_bilatage, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_bilatage", seq_along(.), ""), starts_with("ty_sister_bilatage"))
#variable not collected in bcrpp

#tyrer_data[which (tyrer_data$sisters==1 | tyrer_data$sisters==2),][c(1:50),c("sisters","brcancersis","sister_bc1","sister_bc2","sister_age1","sister_age2",
#                                                       "sister_bilat1","sister_bilat2","sister_bilatage1","sister_bilatage2")]

#still have columns for those with NAs - for all sister and daughter variables created that are NA we can delete the column. e.g., if they have one sister
#they will have NA for sister_bc2, sister_age2, etc. All ending in "2" or more should be deleted for TC dataset.

tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sister_genetic = list({
      n <- sisters_new
      n_bc <-  brcancersis_nomiss
      if (n == 0) {
        #0 if unknown or sisters=0
        rep(NA, sisters_new+1) #Add +1 so variable gets created if 0 sisters
      } 
      else {
        c(rep(0, n_bc), rep(0, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_sister_genetic, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_genetic", seq_along(.), ""), starts_with("ty_sister_genetic"))

############
#OVARIAN DATA
#############

tyrer_data <- tyrer_data%>%
  rowwise() %>%
  mutate(
    ty_sister_ovca = list({
      n <- sisters_new
      n_ov <-  ovcancersis_nomiss
      if (n == 0) {
        rep(NA, sisters_new+1) #Add +1 so variable gets created if 0 sisters
      } 
      else {
        #for sister_bc add a 1 for each sister that was=1
        #add a 0 for each sister that exists that did not have ovca
        c(rep(1, n_ov), rep(0, n - n_ov))
      }
    })
  )%>%
  unnest_wider(ty_sister_ovca, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_ovca", seq_along(.), ""), starts_with("ty_sister_ovca"))




#Not asked in BCRPP  
tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sister_ageov = list({
      n <- sisters_new
      n_ov <-  ovcancersis_nomiss
      if (n == 0) {
        rep(NA, sisters_new+1) #Add +1 so variable gets created if 0 sisters
      } 
      else {
        
        #did not ask for age at ov cancer of sister, so all missing
        c(rep(-99, n_ov), rep(-99, n - n_ov))
      }
    })
  ) %>%
  unnest_wider(ty_sister_ageov, names_sep = "") %>%
  rename_with(~ paste0("ty_sister_ageov", seq_along(.), ""), starts_with("ty_sister_ageov"))  

################
#Daughter info
################

tyrer_data <- tyrer_data%>%
  rowwise() %>%
  mutate(
    ty_daughter_bc = list({
      n <- daughters_new
      n_bc <-  brcancerdau_nomiss
      #if 0 daughters, set as 0
      if (n == 0) {
        rep(NA, daughters_new+1) #add +1 so the variable for BC gets created if 0 daughters
      } 
      else {
        #for daughter_bc add a 1 for each daughter that was=1
        #add a 0 for each daughter that exists that did not have BC (n-n_bc)
        #add NA for NA info on daughters with bc
        #add NA for all the daughters that are not included per individual (i.e., if they have 6 daughters, 6 dummy vars get NA)
        c(rep(1, n_bc), rep(0, n - n_bc))
      }
    })
  )%>%
  unnest_wider(ty_daughter_bc, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_bc", seq_along(.), ""), starts_with("ty_daughter_bc"))

#Not asked in BCRPP
tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_daughter_bilat = list({
      n <- daughters_new
      n_bc <-  brcancerdau_nomiss
      if (n == 0) {
        rep(NA, daughters_new+1)
      } 
      else {
        #0 for all bc unknown
        c(rep(0, n_bc), rep(0, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_daughter_bilat, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_bilat", seq_along(.), ""), starts_with("ty_daughter_bilat"))

#Not asked in BCRPP
tyrer_data<-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_daughter_agebc = list({
      n <- daughters_new
      n_bc <-  brcancerdau_nomiss
      if (n == 0) {
        rep(NA, daughters_new+1)
      } 
      else {
        #for daughter_bc add a 1 for each daughter that was=1
        #add a 0 for each daughter that exists that did not have BC (n-n_bc)
        #add NA for NA info on daughters with bc
        #add NA for all the daughters that are not included per individual (i.e., if they have 6 daughters, 6 dummy vars get NA)
        c(rep(-99, n_bc), rep(-99, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_daughter_agebc, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_agebc", seq_along(.), ""), starts_with("ty_daughter_agebc"))

#variable not collected in bcrpp
tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_daughter_bilatage = list({
      n <- daughters_new
      n_bc <-  brcancerdau_nomiss
      if (n == 0) {
        rep(NA, daughters_new+1)
      } 
      else {
        c(rep(-99, n_bc), rep(-99, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_daughter_bilatage, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_bilatage", seq_along(.), ""), starts_with("ty_daughter_bilatage"))


#variable not collected in bcrpp - 
tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_daughter_genetic = list({
      n <- daughters_new
      n_bc <-  brcancerdau_nomiss
      if (n == 0 ) {
        #0 if unknown or daughters=0
        rep(NA, daughters_new+1)
      } 
      else {
        #for daughter_bc add a 1 for each daughter that was=1
        #add a 0 for each daughter that exists that did not have BC (n-n_bc)
        #add NA for NA info on daughters with bc
        #add NA for all the daughters that are not included per individual (i.e., if they have 6 daughters, 6 dummy vars get NA)
        c(rep(0, n_bc), rep(0, n - n_bc))
      }
    })
  ) %>%
  unnest_wider(ty_daughter_genetic, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_genetic", seq_along(.), ""), starts_with("ty_daughter_genetic"))

tyrer_data <- tyrer_data%>%
  rowwise() %>%
  mutate(
    ty_daughter_ovca = list({
      n <- daughters_new
      n_ov <-  ovcancerdau_nomiss
      if (n == 0) {
        rep(NA,  daughters_new+1)
      } 
      else {
        #for daughter_bc add a 1 for each daughter that was=1
        #add a 0 for each daughter that exists that did not have BC (n-n_bc)
        c(rep(1, n_ov), rep(0, n - n_ov))
      }
    })
  )%>%
  unnest_wider(ty_daughter_ovca, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_ovca", seq_along(.), ""), starts_with("ty_daughter_ovca"))


tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_daughter_ageov = list({
      n <- daughters_new
      n_ov <-  ovcancerdau_nomiss
      if (n == 0) {
        rep(NA,  daughters_new+1)
      } 
      else {
        #did not ask for age at ov cancer of daughter, so all missing
        c(rep(-99, n_ov), rep(-99, n - n_ov))
      }
    })
  ) %>%
  unnest_wider(ty_daughter_ageov, names_sep = "") %>%
  rename_with(~ paste0("ty_daughter_ageov", seq_along(.), ""), starts_with("ty_daughter_ageov"))  

#Additional info


#daughters of sisters
tyrer_data <-tyrer_data %>%
  rowwise() %>%
  mutate(
    ty_sistersdau = list({
      n <- sisters_new
      if (n == 0) {
        rep(NA, sisters_new+1)
      } 
      else {
        #did not ask for age at ov cancer of daughter, so all missing
        c(rep(0, sisters_new))
      }
    })
  ) %>%
  unnest_wider(ty_sistersdau, names_sep = "") %>%
  rename_with(~ paste0("ty_sistersdau", seq_along(.), ""), starts_with("ty_sistersdau"))  


#Now create all TC variables.
#For the dataset- assumption is if they have famhx of Bc=1 listed then it was the mother with BC. Otherwise no assumptions

data_tyrer_assump1 <- tyrer_data%>%
  rowwise() %>%
  mutate(ty_id= gsub("_","",subject_id),
         ty_age=ifelse(!is.na(age), as.integer(floor(age)), -99),
         ty_agemenarche_cont=ifelse(!is.na(agemenarche), floor(agemenarche), -99),
         ty_parity=ifelse(!is.na(parity_new) & parity_new>0, 1, 
                          ifelse(!is.na(parity_new) & parity_new==0, parity_new, 2)),
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
         ty_intended=0, #instructions say put 0 but gui output shows 2 for those with no info on this. questionable.
         ty_timesincehrt=0,
         ty_genetictest=0,
         ty_genetictestfath=0,
         
         #If no data available set to 0
         ty_mother_bc= ifelse(!is.na(brcancermom), brcancermom, 
                              ifelse(is.na(brcancermom) & (!is.na(fhisfstbc)), fhisfstbc, 0)),
         #mother had bilateral bc? unk. for all
         ty_mother_bc_bilat=0,
         #mother ovarian cancer?
         ty_mother_ovarian=ifelse(!is.na(ovcancermom), ovcancermom,
                                  ifelse(is.na(ovcancermom) & (!is.na(fhisfstoc)), fhisfstoc, 0)),
         #age mother had bc -
         #based on distribution of BC cases in US by age SEER Explorer:
         #rates by age 2018-2022: if age <50 the max rate is for 45-50, so set to 45
         # if age >50 the max rate is for 70-75, so set to 70
         # if age not given as <50 or >50 then set to missing
         ty_mother_age_bc=ifelse(ty_mother_bc==1 & fam1grbc50==1 & !is.na(fam1grbc50), 45, 
                                 ifelse(ty_mother_bc==1 & fam1grbc50==0 & !is.na(fam1grbc50), 70,
                                        ifelse(ty_mother_bc==1 & is.na(fam1grbc50), -99, -99))),
         ty_mother_age_bilat_bc=-99, #Not asked in BCRPP
         ty_mother_age_ovar=-99,     #Not asked in BCRPP
         ty_genetictestmoth=0,       #Not asked in BCRPP
         ty_number_sis=sisters_new,
         
         
         ty_patgran=0,           #Not asked in BCRPP
         ty_patgranbilat=0,      #Not asked in BCRPP
         ty_patgranovar=0,       #Not asked in BCRPP
         ty_patgranage=-99,      #Not asked in BCRPP
         ty_patgranagebilat=-99, #Not asked in BCRPP
         ty_patgranageovar=-99,  #Not asked in BCRPP
         ty_patgrangenetic=0,    #Not asked in BCRPP
         ty_matgran=0,           #Not asked in BCRPP
         ty_matgranbilat=0,      #Not asked in BCRPP
         ty_matgranovar=0,       #Not asked in BCRPP
         ty_matgranage=-99,      #Not asked in BCRPP
         ty_matgranagebilat=-99, #Not asked in BCRPP
         ty_matgranageovar=-99,  #Not asked in BCRPP
         ty_matgrangenetic=0,    #Not asked in BCRPP
         
         #for this will have to assume they have 0 paternal and maternal aunts since no info
         ty_pataunts=0,               #Not asked in BCRPP
         #ty_patauntbc=0,              #Not asked in BCRPP
         #ty_patauntbilat=0,           #Not asked in BCRPP
         #ty_patauntovar=0,            #Not asked in BCRPP
         #ty_patauntage=-99,           #Not asked in BCRPP
         #ty_patauntagebilat=-99,      #Not asked in BCRPP  
         #ty_patauntageovar=-99,       #Not asked in BCRPP
         #ty_patauntgenetic=0,         #Not asked in BCRPP
         
         ty_mataunts=0,              #Not asked in BCRPP
         #ty_matauntbc=0,             #Not asked in BCRPP
         #ty_matauntbilat=0,          #Not asked in BCRPP
         #ty_matauntovar=0,           #Not asked in BCRPP
         #ty_matauntage=-99,          #Not asked in BCRPP
         #ty_matauntagebilat=-99,     #Not asked in BCRPP
         #ty_matauntageovar=-99,      #Not asked in BCRPP   
         #ty_matauntgenetic=0,        #Not asked in BCRPP   
         
         #for this will need to assume 0 daughters if no info on daughters
         
         
         ty_daughters=daughters_new,
         #have ty_sisterbc and ty_daughterbc and all those coded.
         #Third stage famhx info
         #No information provided so input zeros
         
         ty_numberofsisplus0=ty_number_sis+0,
         ty_nopathalfsis=0,
         ty_nomathalfsis=0,
         ty_nomataunt=0,
         ty_nopataunt=0,
         ty_nopatuncle=0,
         ty_nomatuncle=0,
         #Last fields- correspond with table 4
         
         ty_fatherbc=ifelse(brcancerdad==1 & !is.na(brcancerdad), 1, 0),
         ty_brotherbc=0,     #Not asked in BCRPP
         ty_density=0,       #Not incorporating for first analysis - most MMD data unavailable
         ty_densityval=-99 ,  
         ty_polyscore= 0)  %>%
  ungroup()#Not yet collected in BCRPP


#Drop columns that aren't needed - this is based on N of sisters or N daughters- so sister_1_bc will be removed if they have 0 sisters.
#order of vars 

#Drop names without ty_

names<-colnames(data_tyrer_assump1)
names_tyrer<-names[grepl("^ty_", names)]

data_tyrer_assump1<-data_tyrer_assump1[,(names(data_tyrer_assump1) %in% names_tyrer)]
#reorder data so the columns are correct
neworder<-c("ty_id","ty_age","ty_agemenarche_cont","ty_parity","ty_agefb","ty_meno","ty_agemeno","ty_heightm","ty_weightkg",
            "ty_hyperplasia","ty_atyp_hyp","ty_lcis","ty_blank","ty_hist_ovarian","ty_age_ovarian","ty_ashkenazi","ty_hrtuse1","ty_hrttype","ty_durhrt",
            "ty_intended","ty_timesincehrt","ty_genetictest","ty_genetictestfath",
            "ty_mother_bc","ty_mother_bc_bilat","ty_mother_ovarian","ty_mother_age_bc","ty_mother_age_bilat_bc","ty_mother_age_ovar","ty_genetictestmoth",
            "ty_number_sis",
            
            "ty_patgran","ty_patgranbilat","ty_patgranovar","ty_patgranage","ty_patgranagebilat","ty_patgranageovar","ty_patgrangenetic",
            "ty_matgran","ty_matgranbilat","ty_matgranovar","ty_matgranage","ty_matgranagebilat","ty_matgranageovar","ty_matgrangenetic",
            "ty_pataunts","ty_mataunts",
            "ty_daughters",
            
            
            "ty_numberofsisplus0", 
           
            "ty_nopathalfsis", "ty_nomathalfsis","ty_nomataunt", "ty_nopataunt","ty_nopatuncle","ty_nomatuncle",
            "ty_fatherbc","ty_brotherbc","ty_density","ty_densityval","ty_polyscore")


data_tyrer_assump1<-data_tyrer_assump1[,neworder]

head(data_tyrer_assump1)
class(data_tyrer_assump1)
#Create a text file that has one line per individual all with different variables based on having NO NA variables. 
#Delete any columns that have NA listed- this will mean columns differ per individual so should do this in a list of dataframes.

#create a new df for each row
row_dataframes_lapply <- lapply(1:nrow(data_tyrer_assump1), function(i) data_tyrer_assump1[i, , drop = FALSE])
head(row_dataframes_lapply)

cleaned_list_of_dataframes <- lapply(row_dataframes_lapply, function(df) {
  df[, colSums(is.na(df)) == 0]
})

#collapse into single line
df_to_single_string<-function(df){
  paste(apply(df, 1, paste, collapse= " "), collapse="")
}


single_line_data<-lapply(cleaned_list_of_dataframes, df_to_single_string)

#NOTE: For external cohorts -- Please replace "N" below with the total number of individuals in the dataset
writeLines(c("v8","N", unlist(single_line_data)), "FILEPATH/data_1to5y_assump1.txt") 

#NOTE: You may have to split up the dataset into chunks based on memory available for running the IBIS risk calculator. 
#The commented lines here provide an example where N total=48935. Generally I was able to get 25,000 ids to work at once.

#writeLines(c("v8","25000", unlist(single_line_data)[1:25000]), "FILEPATH/data_1to5y_assump1_1.txt") 
#writeLines(c("v8","23935", unlist(single_line_data)[25001:48935]), "FILEPATH/data_1to5y_assump1_2.txt") 



###############################################################
# REPEAT THE ABOVE PROCEDURE USING VALDATA2 for 2-6y follow-up
# PLEASE start with valdata2 and repeat all steps above. 
# Write the file as "FILEPATH/data_2to6y_assump1.txt" 
###############################################################









