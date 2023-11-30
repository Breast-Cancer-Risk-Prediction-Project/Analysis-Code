##Model validation comparison in BCRPP

#######################################################
# IF USING BOX --> move on to : "IF NOT USING BOX" otherwise
#######################################################

        ########################################################
        #Step 1. Read in box directory
        ########################################################
        install.packages("tidyverse")
        install.packages("palmerpenguins")
        install.packages(c("boxr", "base", "usethis"))
        
        library(boxr)
        
        ###########################
        # AUTHENTICATION
        ###########################
        # Authentication is the process of syncing RStudio with your Box.com account.
        # If you are not already logged into Box after running this code, you may be
        # asked to log into Box. Identify Box as hard disk in the cloud.
        
        #Replace with your client id and client secret
        
        box_auth(client_id = "clientid",
                 client_secret = "clientsecret") 
        
        # Set the working directory to your Box folder using the folder ID
        # Change dir_id to match the cohort specific directory
        
        box_setwd(dir_id=111111111)
        
        #Create files
        #replace "STUDY" with the name for your cohort
        
        STUDY_core_data =  box_read(file_id = 111111111)
        STUDY_brca_data =  box_read(file_id = 111111111)
        
        STUDY_allcase<-merge(STUDY_core_data, STUDY_brca_data, by="subject_id")
        STUDYcontrol<-setdiff(STUDY_core_data$subject_id,STUDY_brca_data$subject_id)
        STUDY_allcontrol<-STUDY_core_data[which (STUDY_core_data$subject_id %in% STUDYcontrol),]
        
        
        #Match columns in case and control datasets- add columns
        
        names.case<-colnames(STUDY_allcase)
        names.control<-colnames(STUDY_allcontrol)
        #add.names contains a vector of columns from case that you want to add to control dataframe before merging
        add.names <-setdiff(a,b)
        
        for (i in add.names){
          STUDY_allcontrol[,i]<-NA
        }
        STUDY_allcontrol$case<-0


#######################################################
# IF NOT USING BOX:
#######################################################
          
  #######################    
  #Read in your dataset
  #######################
  #These steps assume the datasets are stored separately as core_data.csv and brca_data.csv. If already merged, can skip this step

              STUDY_core_data<-read.csv("filepath/STUDY_core_data.csv")
              STUDY_brca_data<-read.csv("filepath/STUDY_brca_data.csv")
              STUDY_allcase<-merge(STUDY_core_data, STUDY_brca_data, by="subject_id")
              STUDYcontrol<-setdiff(STUDY_core_data$subject_id,STUDY_brca_data$subject_id)
              STUDY_allcontrol<-STUDY_core_data[which (STUDY_core_data$subject_id %in% STUDYcontrol),]
                  
                    #Match columns in case and control datasets- add columns if needed
                    
                    names.case<-colnames(STUDY_allcase)
                    names.control<-colnames(STUDY_allcontrol)
                    #add.names contains a vector of columns from case that you want to add to control dataframe before merging
                    add.names <-setdiff(a,b)
                    
                    for (i in add.names){
                      STUDY_allcontrol[,i]<-NA
                    }
                    STUDY_allcontrol$case<-0
                    
            data1<-rbind(STUDY_allcase, STUDY_allcontrol)

  #######################
  #If core and brca data already merged, read in combined file to start:
  
  data1<-read.csv("filepath/STUDY_combined_data.csv")

##################################
# REQUIRED PROGRAMS
##################################
        if(!require(magrittr)){install.packages("magrittr")}
        library(magrittr)
        
        if(!require(dplyr)){install.packages("dpylr")}
        library(dplyr)
        
        #For plotting calibration
        if(!require(ggplot2)){install.packages("ggplot2")}
        library(ggplot2)

        if(!require(lubridate)){install.packages("lubridate")}
        library(lubridate)

        if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
        BiocManager::install("iCARE", version = 3.18)
        library(iCARE)        


########################
# Formatting variables
########################

#Set dates to NA where needed:
data1<- data1 %>%
  mutate(record_date=ifelse(record_date=="08/08/8000", NA, record_date),
         record_date=ifelse(substr(record_date, 1, 5)=="08/08", sub("^.....","01/01",record_date), record_date),
         record_year=substr(record_date, 7,10),
         #Change dx year to dx year+1 if dx within the first year
         dxdate_primary1new=ifelse(dxdate_primary1==8888, NA, ifelse(dxdate_primary1>record_year, dxdate_primary1, ifelse(dxdate_primary1==record_year, dxdate_primary1+1, NA))),
         dxdate_primary2new=ifelse(dxdate_primary2==8888, NA, ifelse(dxdate_primary2>record_year, dxdate_primary2, ifelse(dxdate_primary2==record_year, dxdate_primary2+1, NA))),
         BBD_year1=ifelse(BBD_year1==8888| BBD_year1==7777, NA, BBD_year1),
         BBD_year2=ifelse(BBD_year2==8888| BBD_year2==7777, NA, BBD_year2),
         BBD_year3=ifelse(BBD_year3==8888| BBD_year3==7777, NA, BBD_year3),
         BBD_year4=ifelse(BBD_year4==8888| BBD_year4==7777, NA, BBD_year4),
         birth_year=ifelse(birth_year==8888, NA, birth_year),
         lastscreen_year=ifelse(lastscreen_year==8888| lastscreen_year==7777| lastscreen_year==888| lastscreen_year==777, NA, lastscreen_year),
         lastfup=ifelse(lastfup==8888, NA, lastfup))

#Change date variables to date format 
#Change record date so that "08/08/YYYY" is stored as 01/01/YYYY
data1 <- data1 %>%
  mutate(record_date=as.Date(record_date, format="%d/%m/%Y"),
         BBD_date1=ymd(BBD_year1, truncated=2L),
         BBD_date2=ymd(BBD_year2, truncated=2L),
         BBD_date3=ymd(BBD_year3, truncated=2L),
         BBD_date4=ymd(BBD_year4, truncated=2L),
         birth_date=ymd(birth_year, truncated=2L),
         lastscreen_date=ymd(lastscreen_year, truncated=2L),
         primary1dx_date=ymd(dxdate_primary1new, truncated=2L),
         primary2dx_date=ymd(dxdate_primary2new, truncated=2L),
         lastfup_date=ymd(lastfup, truncated=2L))

#OUTLIER ASSIGNMENT

data1 <- data1 %>%
  mutate( age_outlier           =ifelse((age<=0| age<777) & age>100, 1, 0),
          height_outlier        =ifelse((height<=0|height<777) & height>200, 1, 0),
          weight_outlier        =ifelse((weight<=0|weight<777)  & weight>150, 1, 0),
          bmi_outlier           =ifelse((bmi<=0   |bmi<777)     & bmi>50, 1, 0),
          waist_outlier         =ifelse((waist<=0|waist<777)    & waist>160, 1, 0),
          hip_outlier           =ifelse((hip<=0|hip<777)        & hip>160, 1, 0),
          wtohip_outlier        =ifelse((whr<=0   |whr<777)     & whr>1.2, 1, 0),
          bmiearly_outlier      =ifelse((bmi_earlyadult<=0|bmi_earlyadult<777) & bmi_earlyadult>50, 1, 0),
          
          alcoholinit_outlier   =ifelse((alcohol_init<=0|alcohol_init<777)     & alcohol_init>70, 1, 0),
          alcoholamt_outlier    =ifelse((alcohol_amt<0|alcohol_amt<777) & alcohol_amt>50, 1, 0),
          #changed this to >50 since the 7 category variable allows up to 45+
          alcoholstop_outlier   =ifelse((alcohol_stop<=0|alcohol_stop<777) & alcohol_stop>100, 1, 0),
          alcoholdur_outlier    =ifelse((alcohol_dur<0|alcohol_dur<777) & alcohol_dur>50, 1, 0),
          smokingdur_outlier    =ifelse((smoking_dur<0|smoking_dur<777) & smoking_dur>50, 1, 0),
          smokingamt_outlier    =ifelse((smoking_amt<0|smoking_amt<777) & smoking_amt>60, 1, 0),
          smokingstop_outlier   =ifelse((smoking_stop<=0|smoking_stop<777) & smoking_stop>100, 1, 0),
          smokinginit_outlier   =ifelse((smoking_init<=0|smoking_init<777) & smoking_init>70, 1, 0),
          agemenarche_outlier   =ifelse((agemenarche<=0|agemenarche<777) & agemenarche>18, 1, 0),
          
          agepreg1_outlier      =ifelse((age_preg1<12|age_preg1<777) & age_preg1>58, 1, 0),
          agepreg2_outlier      =ifelse((age_preg2<12|age_preg2<777) & age_preg2>58, 1, 0),
          agepreg3_outlier      =ifelse((age_preg3<12|age_preg3<777) & age_preg3>58, 1, 0),
          agepreg4_outlier      =ifelse((age_preg4<12|age_preg4<777) & age_preg4>58, 1, 0),
          agepreg5_outlier      =ifelse((age_preg5<12|age_preg5<777) & age_preg5>58, 1, 0),
          agepreg6_outlier      =ifelse((age_preg6<12|age_preg6<777) & age_preg6>58, 1, 0),
          agepreg7_outlier      =ifelse((age_preg7<12|age_preg7<777) & age_preg7>58, 1, 0),
          agepreg8_outlier      =ifelse((age_preg8<12|age_preg8<777) & age_preg8>58, 1, 0),
          agepreg9_outlier      =ifelse((age_preg9<12|age_preg9<777) & age_preg9>58, 1, 0),
          agepreg10_outlier     =ifelse((age_preg10<12|age_preg10<777) & age_preg10>58, 1, 0),
          
          breastfddur1_outlier  =ifelse((breastfeed_dur_b1<0|breastfeed_dur_b1<777) & breastfeed_dur_b1>24, 1, 0),
          breastfddur2_outlier  =ifelse((breastfeed_dur_b2<0|breastfeed_dur_b2<777) & breastfeed_dur_b2>24, 1, 0),
          breastfddur3_outlier  =ifelse((breastfeed_dur_b3<0|breastfeed_dur_b3<777) & breastfeed_dur_b3>24, 1, 0),
          breastfddur4_outlier  =ifelse((breastfeed_dur_b4<0|breastfeed_dur_b4<777) & breastfeed_dur_b4>24, 1, 0),
          breastfddur5_outlier  =ifelse((breastfeed_dur_b5<0|breastfeed_dur_b5<777) & breastfeed_dur_b5>24, 1, 0),
          breastfddur6_outlier  =ifelse((breastfeed_dur_b6<0|breastfeed_dur_b6<777) & breastfeed_dur_b6>24, 1, 0),
          breastfddur7_outlier  =ifelse((breastfeed_dur_b7<0|breastfeed_dur_b7<777) & breastfeed_dur_b7>24, 1, 0),
          breastfddur8_outlier  =ifelse((breastfeed_dur_b8<0|breastfeed_dur_b8<777) & breastfeed_dur_b8>24, 1, 0),
          breastfddur9_outlier  =ifelse((breastfeed_dur_b9<0|breastfeed_dur_b9<777) & breastfeed_dur_b9>24, 1, 0),
          breastfddur10_outlier =ifelse((breastfeed_dur_b10<0|breastfeed_dur_b10<777) & breastfeed_dur_b10>24, 1, 0),
          
          #Set OC use outlier to min age at menarche
          ocusedur_outlier      =ifelse((ocuse_dur<0|ocuse_dur<777) & ocuse_dur>84, 1, 0),
          ocusestart_outlier    =ifelse((ocuse_start<0|ocuse_start<777) & ocuse_start>58, 1, 0),
          ocusestop_outlier     =ifelse((ocuse_stop<=0|ocuse_stop<777) & ocuse_stop>58, 1, 0),
          hrtdur_outlier        =ifelse((hrt_dur<0|hrt_dur<777) & hrt_dur>24, 1, 0),
          hrtepdur_outlier      =ifelse((hrtep_dur<0|hrtep_dur<777) & hrtep_dur>24, 1, 0),
          hrteonlydur_outlier   =ifelse((hrteonly_dur<0|hrteonly_dur<777) & hrteonly_dur>24, 1, 0),
          pamets_outlier        =ifelse((pa_mets<=0|pa_mets<777) & pa_mets>30, 1, 0),
          papct_outlier         =ifelse((pa_pct<0|pa_pct<777) & pa_pct>100, 1, 0),
          
          sizeprimary1_outlier  =ifelse((size_primary1<=0|size_primary1<777) & size_primary1>100, 1, 0),
          sizeprimary2_outlier  =ifelse((size_primary2<=0|size_primary2<777) & size_primary2>100, 1, 0),
          ki67primary1_outlier  =ifelse((ki67_primary1<0|ki67_primary1<777) & ki67_primary1>100, 1, 0),
          ki67primary2_outlier  =ifelse((ki67_primary2<0|ki67_primary2<777) & ki67_primary2>100, 1, 0)
  )


############################################################
# Step 3. Create new variables for risk prediction modeling
############################################################
#Models:
#Gail (BCRAT), CARE, APA (Gail), IBIS, BOADICEA, BWHS
#Create variables based on data dictionary for models
#Treating outliers as NA in models
#Missing variables considered based on model specification

###########################################
# 1. Create variables for dataset 
###########################################

#Fix so parous and age at first birth align
#Fix so alcohol current gm and status align
#Create categorical variables based on different risk models
#Fix Biopsy and Atypia so if biopsy is 0 then Atypia has to be missing (req. in Gail model)


data2<- data1 %>% 
  mutate( ethnicity=ifelse(ethnicity==888, NA, ethnicity),
          race=ifelse(race !=888, race, NA),
          age=ifelse(age<777 & age_outlier !=1, age, NA),
          agemenarche_cont=ifelse(agemenarche<777 & agemenarche_outlier !=1, agemenarche, NA),
          agemenarche_cat3=ifelse(agemenarche_cont>0 & agemenarche_cont<12, 1, 
                                  ifelse(agemenarche_cont >=12 & agemenarche_cont<14, 2, 
                                         ifelse(agemenarche_cont >=14, 3, NA))), 
          agemenarche_cat7=ifelse(agemenarche_cont>0 & agemenarche_cont<10.5, 1, 
                                  ifelse(agemenarche_cont>=10.5 & agemenarche_cont<11.5, 2, 
                                         ifelse(agemenarche_cont>=11.5 & agemenarche_cont<12.5, 3, 
                                                ifelse(agemenarche_cont>=12.5& agemenarche_cont<13.5, 4, 
                                                       ifelse(agemenarche_cont>=13.5&agemenarche_cont<14.5, 5, 
                                                              ifelse(agemenarche_cont>=14.5&agemenarche_cont<15.5, 6, 
                                                                     ifelse(agemenarche_cont>=15.5, 7, NA))))))),
          agemenarche_cat2=ifelse(agemenarche_cont>0 & agemenarche_cont<14, 1, 
                                  ifelse(agemenarche_cont>=14, 2, NA)),
          age_ge50=ifelse(age<50 & !is.na(age), 0, 
                          ifelse(age>=50, 1, NA)),
          age_cat1=ifelse(age>=45 & age<50, 1, 
                          ifelse(age>=50 & age<60,2,NA)),
          age_cat2=ifelse(age>=45 & age<50, 1, 
                          ifelse(age>=50&age<55, 2, 
                                 ifelse(age>=55&age<60, 3, 
                                        ifelse(age>=60 & age<65, 4, 
                                               ifelse(age>=65 & age<75, 5, NA))))),
          meno_status=ifelse(meno_status==888, NA, meno_status),
          meno_stat=ifelse(meno_status==2 | meno_status==3, 1, ifelse(meno_status==1, 2, NA)),
          age_menostat=ifelse(meno_stat==1 & age_cat1==1, 1, 
                              ifelse(meno_stat==1 & age_cat1==2, 2, 
                                     ifelse(meno_stat==2 & age_cat2==1, 3, 
                                            ifelse(meno_stat==2 & age_cat2==2, 4, 
                                                   ifelse(meno_stat==2 &age_cat2==3, 5, 
                                                          ifelse(meno_stat==2 &age_cat2==4, 6,NA)))))),
          parous=ifelse(parous==888, NA, parous),
          ageflb=ifelse(age_preg1<777 & age_preg1>=agemenarche_cont & !is.na(agemenarche_cont) & age_preg1<=age & !is.na(age) & agepreg1_outlier !=1 , age_preg1, NA),
          #recode parous if parous and ageflb don't match - change parous to 1 if there is a valid age at first birth that fits in range above
          parous_new=ifelse(!is.na(ageflb), 1, parous),
          parity_new=ifelse(parous_new==0, 0,
                            ifelse(parous_new==1 & parity>0 & parity<777, parity, 
                                   ifelse(parous_new==1 & (parity==0 | parity==888|parity==777), 1, NA))),
          parity_cat=ifelse(parity_new<3, parity_new, 3),
          ageflb_cat5=ifelse(parous_new==0, 1, 
                             ifelse(ageflb>0 & ageflb<20, 2, 
                                    ifelse(ageflb>=20 & ageflb<25, 3, 
                                           ifelse(ageflb>=25 & ageflb<30, 4, 
                                                  ifelse(ageflb>=30, 5, NA))))),
          ageflb_cat4=ifelse(parous_new==0, 1, 
                             ifelse(ageflb>0 &ageflb<25, 2, 
                                    ifelse(ageflb>=25&ageflb<30, 3, 
                                           ifelse(ageflb>=30, 4, NA)))),
          ageflb_parous=ifelse(ageflb>0 & ageflb<20, 1, 
                               ifelse(ageflb>=20 & ageflb<25, 2, 
                                      ifelse(ageflb>=25&ageflb<30, 3, 
                                             ifelse(ageflb>=30, 4, NA)))),
          famhx_first=ifelse(fhx_fdr_brca==0, 0, 
                             ifelse(fhx_fdr_brca==1, 1, NA)),
          ageflb_famhx=ifelse(ageflb_cat5==1&famhx_first==0, 1, 
                              ifelse(ageflb_cat5==1&famhx_first==1, 2, 
                                     ifelse(ageflb_cat5==2&famhx_first==0, 3, 
                                            ifelse(ageflb_cat5==2&famhx_first==1, 4, 
                                                   ifelse(ageflb_cat5==3&famhx_first==0, 5, 
                                                          ifelse(ageflb_cat5==3&famhx_first==1, 6, 
                                                                 ifelse(ageflb_cat5==4&famhx_first==0, 7, 
                                                                        ifelse(ageflb_cat5==4&famhx_first==1, 8, 
                                                                               ifelse(ageflb_cat5==5&famhx_first==0, 9, 
                                                                                      ifelse(ageflb_cat5==5&famhx_first==1, 10, NA)))))))))),
          biopsy_num=ifelse(Biopsies_number !=888, Biopsies_number, 
                            ifelse(Biopsies_number==888 & Biopsies_yesno==1, 1, 
                                   ifelse(Biopsies_number==888 & Biopsies_yesno==0,0, NA))),
          atyp_hyp=ifelse(biopsy_num>0 & (BBD_type1==3 | BBD_type2==3 | BBD_type3==3 | BBD_type4)==3, 1, 
                          ifelse(biopsy_num>0 & (BBD_type1==1|BBD_type1==2|BBD_type2==1|BBD_type2==2|BBD_type3==1|BBD_type3==2), 0, NA)),
          biop_noatyp=ifelse(biopsy_num>0 & atyp_hyp==0, 1, 
                             ifelse(biopsy_num>0 & atyp_hyp==1, 0, NA)),
          bbd=ifelse(Biopsies_yesno !=888, Biopsies_yesno, NA),
          biopsy_cat3=ifelse(Biopsies_yesno==0 | BBD_type1==1, 1, 
                             ifelse(Biopsies_yesno==1 & BBD_type1 !=2, 2, 
                                    ifelse(BBD_type1==2, 3, NA))),
          biopsy_age50=ifelse(biopsy_num==0 & age_ge50==0, 1, 
                              ifelse(biopsy_num==1&age_ge50==0, 2, 
                                     ifelse(biopsy_num>1 & age_ge50==0, 3, 
                                            ifelse(biopsy_num==0 & age_ge50==1, 4, 
                                                   ifelse(biopsy_num==1 & age_ge50==1, 5, 
                                                          ifelse(biopsy_num>1&age_ge50==1, 6, NA)))))),
          postmeno_dur_hold=ifelse (meno_age!=777 & meno_age !=888 & meno_stat==2, age-meno_age, 
                                    ifelse(meno_stat==1, 0, NA)),
          postmeno_dur=ifelse(postmeno_dur_hold<0, 0, postmeno_dur_hold),
          hrt_use=ifelse(hrtuse==0|meno_stat==1, 1, 
                         ifelse(hrtuse==1 & meno_stat==2, 2, 
                                ifelse(hrtuse==2 & meno_stat==2, 3, 
                                       ifelse(hrtuse==3 & meno_stat==2, 4, NA)))),
          hrt_use_icare=ifelse(hrt_use<4, hrt_use, 3),
          bmi_earlyadult=ifelse(bmi_earlyadult<777 & bmiearly_outlier !=1, bmi_earlyadult, NA),  
          bmi_cont=ifelse(bmi<777 & bmi_outlier !=1, bmi, NA),
          bmi_cat3=ifelse(bmi_cont>0 & bmi_cont<25, 1,
                          ifelse(bmi_cont>=25 & bmi_cont<30, 2, 
                                 ifelse(bmi_cont>=30, 3, NA))),
          bmi_cat4=ifelse(bmi_cont>0 & bmi_cont<18.5, 1, 
                          ifelse(bmi_cont>=18.5&bmi_cont<25, 2, 
                                 ifelse(bmi_cont>=25&bmi_cont<30, 3, 
                                        ifelse(bmi_cont>=30, 4, NA)))),
          bmi_cat5=ifelse(bmi_cont>0 & bmi_cont<21, 1, 
                          ifelse(bmi_cont>=21&bmi_cont<23, 2, 
                                 ifelse(bmi_cont>=23&bmi_cont<35, 3, 
                                        ifelse(bmi_cont>=25 & bmi_cont<27, 4, 
                                               ifelse(bmi_cont>=27, 5, NA))))),
          bmi_cat2=ifelse(bmi_cont>0 & bmi_cont<30, 1, 
                          ifelse(bmi_cont>=30,2, NA )),
          pmh_type_cur=ifelse(hrtuse==0|hrtuse==2, 1, 
                              ifelse(hrtuse_eonly==1, 3, 
                                     ifelse(hrtuse_ep==1, 4,
                                            ifelse (hrtuse==1, 2, NA)))),
          pmh_type_icare=ifelse(pmh_type_cur==3, 0, ifelse(pmh_type_cur==4, 1, NA)), 
          pmh_type_any=ifelse(hrtuse==0, 1, 
                              ifelse((hrtuse==2|hrtuse==3) & hrtuse_eonly!=2 & hrtuse_eonly!=3 & hrtuse_ep!=2 & hrtuse_ep!=3, 2, 
                                     ifelse((hrtuse==2 | hrtuse==3) & (hrtuse_eonly==2 | hrtuse_eonly==3), 3, 
                                            ifelse((hrtuse==2 | hrtuse==3) & (hrtuse_ep==2|hrtuse_ep==3), 4, 
                                                   ifelse(pmh_type_cur==2, 5, 
                                                          ifelse(pmh_type_cur==3, 6, 
                                                                 ifelse(pmh_type_cur==4, 7, NA))))))),
          height_cont=ifelse(height<777 & height_outlier !=1, height, NA),
          height_in=height_cont/2.54,
          height_cat3=ifelse(height_cont>0 & height_cont <160, 1,
                             ifelse(height_cont >=160 & height_cont<170, 2, 
                                    ifelse(height_cont>=170, 3, NA))),
          #reclassify alcohol status if "unkn. current or former and has alcohol_amt list for current use, change to current
          alcohol_status=ifelse(alcohol_amt>0 & alcohol_amt<777 & alcohol_status==3, 1, alcohol_status),
          #Only assign an amount if they listed current drinker or unknown current/never, otherwise 0 if not current drinker, NA if unk. 
          alcgm=ifelse(alcohol_amt==777 | alcohol_status==2 | alcohol_status==4, 0,
                       ifelse((alcohol_status==1|alcohol_status==3) & alcohol_amt<777 & alcoholamt_outlier !=1, alcohol_amt, NA)),
          alcgm_cat4=ifelse(alcohol_status==2|alcohol_status==4, 1, 
                            ifelse(alcgm>0 & alcgm<11, 2,
                                   ifelse(alcgm>=11 & alcgm<22, 3,
                                          ifelse(alcgm>=22, 4, NA)))),
          alcgm_cat7=ifelse(alcgm_cat4==1, 1, 
                            ifelse(alcgm>0& alcgm<5, 2,
                                   ifelse(alcgm>=5&alcgm<15,3,
                                          ifelse(alcgm>=15 & alcgm<25, 4,
                                                 ifelse(alcgm>=25&alcgm<35, 5,
                                                        ifelse(alcgm>=25 & alcgm<45, 6,
                                                               ifelse(alcgm>=45, 7, NA))))))),
          ocuse=ifelse(ocuse_ever==0,1,
                       ifelse(ocuse_ever==1 & ocuse_current==888, 2,
                              ifelse(ocuse_ever==1 & ocuse_current==0, 3,
                                     ifelse(ocuse_current==1, 4, NA)))),
          meno_age=ifelse(meno_age<777, meno_age, NA),
          meno_age_cat5_icare=ifelse(is.na(meno_age), NA, ifelse(meno_age>=50 & meno_age<55, 1, ifelse(meno_age<40, 2, ifelse(meno_age>=40 & meno_age<45, 3, ifelse(meno_age>=45 & meno_age<50, 4, 5))))),
          bca_case=ifelse(case==0, 0, 
                          ifelse(case==1, 1, NA)),
          bca_invasive=ifelse(bca_case==1 & invasive_primary1==1, 1, 
                              ifelse(bca_case==1 & invasive_primary1==2, 0, NA)),
          censor_year=lastfup,
          #change missings to NA for alcohol and smoking status
          alcohol_status=ifelse(alcohol_status==888, NA, alcohol_status),
          smoking_status=ifelse(smoking_status==888, NA, smoking_status),
          #if censor_date equal to record year then change to be record year + 1
          censor_year=ifelse(censor_year==8888, NA, ifelse(censor_year==record_year, censor_year+1, censor_year)),
          censor_date=ymd(censor_year, truncated=2L),
          timediff_control=as.numeric(censor_date-record_date),
          timediff_control_yr=as.numeric(timediff_control/365),
          timediff_case=as.numeric(primary1dx_date-record_date),
          timediff_case_yr=as.numeric(timediff_case/365),
          case_5y=ifelse(case==1 & timediff_case_yr<=5 & !is.na(timediff_case_yr), 1, 0),
          case_10y=ifelse(case==1 & timediff_case_yr<=10 & !is.na(timediff_case_yr), 1, 0))


#######################
# Add labels- 
#######################

install.packages("expss")
library(expss)
data2=apply_labels(data2,
                   ethnicity="Ethnicity",
                   ethnicity=c("Nonhispanic"=0, "Hispanic"=1, "NA"=NA ),
                   bmi="BMI, kg/m2",
                   bmi_earlyadult="BMI in early adulthood, kg/m2",
                   agemenarche_cont="Age at menarche, years",
                   age="Age,years",
                   agemenarche_cat3="Age at menarche, Gail",
                   agemenarche_cat3=c("NA"=NA,"<12"=1, "12-<14"=2, ">=14"=3),
                   agemenarche_cat7="Age at menarche, iCARE BOAD",
                   agemenarche_cat7=c("NA"=NA, "<10.5"=1, "10.5-<11.5"=2, "11.5-<12.5"=3, "12.5-<13.5"=4, "13.5-<14.5"=5, "14.5-<15.5"=6, "15.5-16"=7),
                   agemenarche_cat2="Age at menarche, BWHS",
                   agemenarche_cat2=c("NA"=NA, "<14"=1, "14-16"=2),
                   biopsy_num="Number of biopsies",
                   Biopsies_yesno="Had biopsy",
                   age_ge50="Age >=50",
                   age_ge50=c("NA"=NA, "<50"=0, ">=50-100"=1),
                   biopsy_age50="Biopsy and age combined",
                   biopsy_age50=c("NA"=NA, "0 biopsies, <50"=1, "1 biopsy, <50"=2,"2 biopsies, <50"=3,"0 biopises, >=50"=4,"1 biopsy, >=50"=5, "2 biopsies, >=50"),
                   ageflb="Age at first pregnancy",
                   ageflb_cat5="Age at first birth, Gail, TC",
                   ageflb_cat5=c("NA"=NA, "Non-parous"=1, "<20"=2, "20-<25"=3, "25-<30"=4, "30-58"=5),
                   ageflb_cat4="Age at first birth, Rosner",
                   ageflb_cat4=c("NA"=NA, "Non-parous"=1, "<25"=2, "25-<30"=3, "30-58"=4),
                   ageflb_parous="Age at first birth if parous, BOAD, iCARE",
                   ageflb_parous=c("NA"=NA, "<20"=1, "20-<25"=2, "25-<30"=3, ">=30"=4),
                   parity_cat="Parity, categorical",
                   parity_cat=c("NA"=NA, "0"=0, "1"=1, "2"=2, ">=3"=3),
                   famhx_first="First-degree family history of breast cancer",
                   famhx_first=c("NA"=NA, "No"=0, "Yes"=1),
                   ageflb_famhx="Age at first birth and famhx combined",
                   ageflb_famhx=c("NA"=NA, "Nonparous, no famhx"=1, "Nonparous, yes famhx"=2, "<20, no famhx"=3, "<20, yes famhx"=4,"20-<25, no famhx"=5,"20-<25, yes famhx"=6, "25-<30 no famhx"=7, "25-<30, yes famhx"=8, "30+, no famhx"=9, "30+, yes famhx"=10),
                   atyp_hyp="Atypical hyperplasia",
                   atyp_hyp=c("NA"=NA, "No"=0, "Yes"=1),
                   biop_noatyp="Had biopsy, but not atypia",
                   biop_noatyp=c("NA"=NA, "No"=0, "Yes"=1),
                   postmeno_dur="Duration of postmenopause in years",
                   age_cat1="Age group, 45-60",
                   age_cat1=c("NA"=NA,"45-<50"=1,"50-<60"=2),
                   age_cat2="Age group, 45-75",
                   age_cat2=c("NA"=NA,"45-<50"=1, "50-<55"=2, "55-<60"=3, "60-<65"=4, "65-<75"=5),
                   meno_stat="Menopausal status pre post",
                   meno_stat=c("NA"=NA,"Premenopausal"=1,"Postmenopausal"=2),
                   meno_age="Age at menopause, years",
                   meno_age_cat5_icare="Age at menopause, categories",
                   meno_age_cat5_icare=c("NA"=NA, "50-55"=1,"<40"=2,"40-45"=3,"45-50"=4,">=55"=5),
                   age_menostat="Age and menopausal status combined",
                   age_menostat=c("Pre, 45-<50"=1, "Pre, 50<60"=2,"Post, 45-<50"=3, "Post, 50-<55"=4, "Post, 55-<60"=5, "Post, 60-<65"=6),
                   bmi_cat3="BMI categorical, iCARE, ge50",
                   bmi_cat3=c("NA"=NA, "<25"=1, "25-<30"=2, ">30-50"=3),
                   bmi_cat4="BMI categorical, iCARE, BOAD",
                   bmi_cat4=c("NA"=NA, "<18.5"=1, "18.5-<25"=2, "25-<30"=3,"30-50"=4),
                   bmi_cat5="BMI categorical, TC",
                   bmi_cat5=c("NA"=NA, "<21"=1, "21-<23"=2, "23-<25"=3, "25-<27"=4, "27-50"=5),
                   bmi_cat2="BMI under over 30, BWHS",
                   bmi_cat2=c("NA"=NA, "<30"=1, "30-50"=2),
                   bbd="BBD number, rosner",
                   biopsy_cat3="Biopsy and BBD type, TC",
                   biopsy_cat3=c("NA"=NA, "No biopsy or no proliferative disease"=1, "Prior biopsy, result unknown"=2,"Hyperplasia (not atypia)"),
                   hrt_use="HRT use",
                   hrt_use=c("NA"=NA, "Never use"=1, "Current use"=2, "Former use"=3, "Ever use (unk)"=4),
                   hrt_use_icare=c("NA"=NA, "never"=1, "current"=2, "former"=3),
                   pmh_type_cur="Type current PMH, Rosner",
                   pmh_type_cur=c("NA"=NA, "None"=1, "Unknown type"=2,"Estrogen only"=3, "Combined"=4),
                   pmh_type_any="Ever use and type PMH",
                   pmh_type_any=c("NA"=NA, "Never"=1, "Past, unk type"=2, "Past, E only"=3, "Past, Combined"=4, "Current, unk type"=5, "Current, E only"=6, "Current, combined"=7),
                   pmh_type_icare=c("NA"=NA, "Estrogen only"=0, "Combined"=1),
                   height="Height, cm",
                   height_in="Height, in",
                   height_cat3="Height categories, TC",
                   height_cat3=c("NA"=NA, "0-<160cm"=1, "160-<170cm"=2, "170-200cm"=3),
                   alcgm_cat4="Alcohol grams, categorical",
                   alcgm_cat4=c("NA"=NA,"None"=1, "<11g"=2, "11-<22g"=3, "22-<50g"=4),
                   alcgm_cat7="Alcohol grams, 7 categories",
                   alcgm_cat7=c("NA"=NA, "0"=1, "<5"=2, "5-<15"=3, "15-<25"=4, "25-<35"=5, "35-<45"=6, "45-50"=7),
                   ocuse="OC use",
                   ocuse=c("NA"=NA,"Never"=1,"Ever, unk. past or current"=2, "Past"=3, "Current"=4),
                   dxdate_primary1="Date of first primary dx",
                   bca_case="Case status",
                   bca_case=c("NA"=NA, "Non-case"=0, "Case"=1),
                   censor_date="Date of censoring",
                   bca_invasive="Invasive or in situ",
                   bca_invasive=c("NA"=NA, "Invasive"=1, "In situ"=0),
                   alcohol_status="Alcohol status",
                   alcohol_status=c("Current"=1,"Former"=2,"Ever, unk. if current or former"=3, "Never"=4,"Missing"=NA),
                   smoking_status="Smoking status",
                   smoking_status=c("Current"=1,"Former"=2,"Ever, unk. if current or former"=3, "Never"=4,"Missing"=NA),
                   meno_status="Menopausal status all",
                   meno_status=c("Postmenopausal"=1, "Premenopausal"=2, "Perimenopausal/other"=3, "Missing"=NA),
                   ethnicity=c("NA"=NA, "Hispanic"=1, "Nonhispanic"=0),
                   race=c("White"=1, "Black"=2, "Asian"=3,"Native Hawaiian/Pacific Islander"=4, "American Indian, Alaska Native"=5, "Other/multiracial"=6, "Missing"=NA),
                   parous_new=c("Nonparous"=0,"Parous"=1, "Missing"=NA),
                   bbd=c("No BBD"=0, "Yes BBD"=1, "NA"=NA))


###############################
# Follow-up information
##############################

data2 <- data2 %>%
  mutate(dxtofup= lastfup_date-primary1dx_date,
         dxtofup= as.integer(dxtofup),
         dxtofup_yr=dxtofup/365.25,
         dxage=primary1dx_date-birth_date,
         dxage= as.integer(dxage),
         dxage_yr=dxage/365.25)
#summary(data2$dxage_yr[which (data2$case==1)])        

data2 <-data2 %>%
  mutate(size_primary1_new=ifelse(sizeprimary1_outlier==1 | size_primary1==888, NA, size_primary1),
         ki67_primary1_new=ifelse(ki67primary1_outlier==1 | ki67_primary1==888, NA, ki67_primary1),
         ki67cat_primary1=ifelse(ki67cat_primary1==888 | ki67cat_primary1==777, NA, ki67cat_primary1),
         sizecat_primary1=ifelse(sizecat_primary1==888 | sizecat_primary1==777, NA, sizecat_primary1),
         her2_primary1=ifelse(her2_primary1==888, NA, her2_primary1),
         er_primary1=ifelse(er_primary1==888 | er_primary1==777, NA, er_primary1),
         pr_primary1=ifelse(pr_primary1==888 | pr_primary1==777, NA, pr_primary1),
         grade_primary1=ifelse(grade_primary1==888 | grade_primary1==777, NA, grade_primary1),
         stage_primary1=ifelse(stage_primary1==888 | stage_primary1==777, NA, stage_primary1))


###########################################
# PRINT DESCRIPTIVE STATS (if desired)- 
##########################################

#This section will just print out descriptive stats for variables of interest to double check that the coding looks correct.

data2_post<-data2[which (data2$meno_stat==2),]

data2_continuous<-data2[c("age",
                          "agemenarche_cont",
                          "ageflb",
                          "parity_new",
                          "postmeno_dur",
                          "meno_age",
                          "bmi_cont",
                          "height_in",
                          "bmi_earlyadult","case")]

data2_continuouspost<-data2_post[c("age",
                                   "agemenarche_cont",
                                   "ageflb",
                                   "parity_new",
                                   "postmeno_dur",
                                   "meno_age",
                                   "bmi_cont",
                                   "height_in",
                                   "bmi_earlyadult","case")]
#colnames(data2)
data2_categories<-data2[c("smoking_status",
                          "alcohol_status",
                          "biopsy_num",
                          "ethnicity",
                          "sex",
                          "race",
                          "age_cat1",
                          "age_cat2",
                          "meno_status",
                          "bmi_cat2",
                          "bmi_cat4",
                          "bmi_cat5",
                          "agemenarche_cat3",
                          "agemenarche_cat7",
                          "agemenarche_cat2",
                          "age_ge50",
                          "biopsy_age50",
                          "parous_new",
                          "ageflb_cat5",
                          "ageflb_cat4",
                          "ageflb_parous",
                          "famhx_first",
                          "ageflb_famhx",
                          "atyp_hyp",
                          "biop_noatyp",
                          "bbd",
                          "biopsy_cat3",
                          "pmh_type_cur",
                          "pmh_type_any",
                          "height_cat3",
                          "alcgm_cat4",
                          "alcgm_cat7",
                          "ocuse","case")]

data2_categoriespost<-data2_post[c("smoking_status",
                                   "alcohol_status",
                                   "biopsy_num",
                                   "ethnicity",
                                   "sex",
                                   "race",
                                   "age_cat1",
                                   "age_cat2",
                                   "meno_status",
                                   "bmi_cat2",
                                   "bmi_cat4",
                                   "bmi_cat5",
                                   "agemenarche_cat3",
                                   "agemenarche_cat7",
                                   "agemenarche_cat2",
                                   "age_ge50",
                                   "biopsy_age50",
                                   "parous_new",
                                   "ageflb_cat5",
                                   "ageflb_cat4",
                                   "ageflb_parous",
                                   "famhx_first",
                                   "ageflb_famhx",
                                   "atyp_hyp",
                                   "biop_noatyp",
                                   "bbd",
                                   "biopsy_cat3",
                                   "pmh_type_cur",
                                   "pmh_type_any",
                                   "height_cat3",
                                   "alcgm_cat4",
                                   "alcgm_cat7",
                                   "ocuse","case")]

continuous_summary<-cbind(min=apply(data2_continuous,2,min,na.rm=TRUE), median=apply(data2_continuous,2,median,na.rm=TRUE),max=apply(data2_continuous,2,max,na.rm=TRUE),IQR=apply(data2_continuous,2,IQR,na.rm=TRUE), mean=apply(data2_continuous,2,mean,na.rm=TRUE),sd=apply(data2_continuous,2,sd,na.rm=TRUE))
continuous_summary_post<-cbind(min=apply(data2_continuouspost,2,min,na.rm=TRUE), median=apply(data2_continuouspost,2,median,na.rm=TRUE),max=apply(data2_continuouspost,2,max,na.rm=TRUE),IQR=apply(data2_continuouspost,2,IQR,na.rm=TRUE), mean=apply(data2_continuouspost,2,mean,na.rm=TRUE),sd=apply(data2_continuouspost, 2,sd,na.rm=TRUE))

#Replace filepaths if you want to write as tables.

write.table(continuous_summary, file="filepath/continuousvar_summary_STUDY.txt", sep=',')
write.table(continuous_summary_post, file="filepath/continuousvar_summary_post_STUDY.txt", sep=',')

#Categorical Variables

tblFun <- function(x){
  tbl <- table(x, useNA="always")
  res <- cbind(tbl,round(prop.table(tbl)*100,2))
  colnames(res) <- c('Count','Percentage')
  res
}

summariescat<-lapply(data2_categories[1:34],tblFun)
capture.output(summariescat, file="filepath/categorical_summary_STUDY.txt",sep=',')

summariescat_post<-lapply(data2_categoriespost[1:34],tblFun)
capture.output(summariescat_post, file="filepath/categorical_summary_post_STUDY.txt", sep=',')

    
    ##############################################
    # CASE DESCRIPTIVE STATS
    ##############################################

      data2_case<-data2[which (data2$case==1),]
  
      data2_case_cont<-data2_case[c("dxage_yr", "dxtofup_yr", "size_primary1_new", "ki67_primary1_new")]
      data2_case_cat<-data2_case[c("stage_primary1", "grade_primary1", "sizecat_primary1","er_primary1","pr_primary1","her2_primary1","ki67cat_primary1")]
  
      case_summary<-cbind(min=round(apply(data2_case_cont,2,min,na.rm=TRUE), digits=3), median=apply(data2_case_cont,2,median,na.rm=TRUE),max=apply(data2_case_cont,2,max,na.rm=TRUE),IQR=apply(data2_case_cont,2,IQR,na.rm=TRUE), mean=apply(data2_case_cont,2,mean,na.rm=TRUE),sd=apply(data2_case_cont,2,sd,na.rm=TRUE))
      write.table(case_summary, file="filepath/case_summary1_STUDY.txt", sep=',')


      case_summary2<-lapply(data2_case_cat[1:7],tblFun)
      capture.output(case_summary2, file="filepath/case_summary2_STUDY.txt", sep=',')


#####################################################################################
# GAIL Model - Calculate Relative and Absolute risk for calibration
#####################################################################################

#In Gail model, all missing variables are set to the most common

#Beta coefficients
#Load from R package available online
install.packages("BCRA")
library(BCRA)
#has an existing code to "estimaterelative risks"

#Create dataset, variable names & assignments for Gail model 
data_gail <-data2[,c("id","age","race","ethnicity","agemenarche_cont","ageflb","famhx_first","biopsy_num","atyp_hyp","parous_new")]

#Reassign values based on Gail parameters

#Missing data assignments
#The standards used in the online Gail model were used for variables with missing data. 
#Women with unknown age at menarche were assigned menarche at 14 years of age or older. 
#Women with unknown age at first live birth were classified as giving birth before age 20 years. 
#Women with missing family history were classified as having no family history. 
#To verify agreement between the code we used and the online tool, we randomly selected 10 patients 
#and compared the 5-year and lifetime risk estimates obtained from the code given us to those from the online 
#risk assessment tool. All of the estimates were in complete agreement.

table(data_gail$race, data_gail$ethnicity)
#For our purposes, assign all Asian under "Chinese" for abs.risk calibration
#Consider all hispanics american born
#1=Wh White 1983-87 SEER rates (rates used in NCI BCRAT)
#2=AA African-American
#3=HU Hispanic-American (US born) 1995-04
#4=NA Other (Native American and unknown race)
#5=HF Hispanic-American (Foreign born) 1995-04
#6=Ch Chinese-American
#7=Ja Japanese-American
#8=Fi Filipino-American
#9=Hw Hawaiian-American/*
#10=oP Other Pacific Islander
#11=oA Other Asian


data_gail2 <- data_gail %>%
  mutate(gailrace=ifelse(race==1 & ethnicity==1, 3, ifelse(race==1 & ethnicity !=1 ,1, ifelse(race==2, 2, ifelse(race==5, 4, ifelse(race==3, 6, ifelse(race==4, 9, 4)))))),
         gailrace=replace(gailrace, is.na(gailrace), 4),
         gailagemen=ifelse(is.na(agemenarche_cont), 99, agemenarche_cont),
         gailageflb=replace(ageflb, is.na(ageflb) & parous_new==0, 98),
         gailageflb=replace(gailageflb, is.na(gailageflb), 99),
         gailnrels=ifelse(is.na(famhx_first), 99, famhx_first),
         gailbiop=ifelse(is.na(biopsy_num), 99, biopsy_num),
         gailatyp=ifelse(is.na(atyp_hyp), 99, atyp_hyp),
         gailageend=age+5)


data_gail_final<-data_gail2[,c("id","age","gailageend","gailrace","gailagemen","gailageflb","gailnrels","gailbiop","gailatyp")]
#Change names
names(data_gail_final)<-c("id","T1","T2","Race","AgeMen","Age1st","N_Rels","N_Biop","HypPlas")

#Data check! 
#Make sure the variables are recoded correctly and check errors
recodechecktable<-recode.check(data_gail_final, Raw_Ind=1)
test<-recodechecktable[which (recodechecktable$Error_Ind>0),]
test
#If any observations show up on this print out then these have errors that need to be fixed- in the assignment of raw gail vars
#No errors

#Estimate rel risk
relrisk_gail<-relative.risk(data_gail_final, Raw_Ind=1)
data_gail_final2<-cbind(data_gail_final, relrisk_gail)
#RR_Star1=relative risk <50
#RR_Star2=RR for woman >=50

#Assign rel risk based on age <50 or >50 to dataset
data_gail_final2$Relrisk_gail<-ifelse(data_gail_final2$T1>=50, data_gail_final2$RR_Star2, data_gail_final2$RR_Star1)

#Now apply Absolute risk
absrisk_gail<-absolute.risk(data_gail_final2, Raw_Ind=1, Avg_White=0)

#Assign absrisk in dataset
data_gail_final2$absrisk_gail<-absrisk_gail

#Data check!
#Make sure mean for Error_Ind=0, otherwise there is an error in the file
#Note, abs risk is given as a % for 5 years
check.summary(data_gail_final, Raw_Ind=1, Avg_White=0)

# > check.summary(data_gail_final, Raw_Ind=1, Avg_White=0)
# Variable                                       Label              Mean            StdDev      N NMiss
# 1 Error_Ind        If mean not 0, implies ERROR in file                 0                 0 314062     0
# 2   AbsRisk Abs risk(%) of BrCa in age interval [T1,T2) 0.797381255458494 0.642943913035021 314062     0
# 3  RR_Star1                     Relative risk age lt 50  1.69609362957351 0.532036284658112 314062     0
# 4  RR_Star2                     Relative risk age ge 50  1.67885275489836 0.494068534492978 314062     0


#For comparison, classify individuals into relative risk categories
#<15%, 15-20%, >20%
data_gail_final2$relriskcat<-ifelse(data_gail_final2$Relrisk_gail<1.15, 1, ifelse(data_gail_final2$Relrisk_gail>=1.15 & data_gail_final2$Relrisk_gail<1.54, 2, 3))
data_gail_final2$absrisk_gailnew<-data_gail_final2$absrisk_gail/100


################
# Check distribution of absolute and relative risks
#Plot distribution of relative risk in Gail model
# 
# ggplot(data_gail_final2, aes(x=Relrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Relative Risk",y="Frequency") + ggtitle("Relative Risk Distribution for Gail Model")
# ggsave("filepath/Gail_relrisk.png", width=4, height=4)
# 
# ggplot(data_gail_final2, aes(x=absrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Absolute Risk",y="Frequency") + ggtitle("Relative Risk Distribution for Gail Model")
# ggsave("filepath/Gail_absrisk.png", width=4, height=4)
# 
# summary(data_gail_final2$Relrisk_gail)

# #Using abs. risk / 100
# summary(data_gail_final2$absrisk_gailnew)

#########################

#Subset dataset by race as determined in Gail model

data_gail_final_white<-data_gail_final2[which (data_gail_final2$Race==1),] #Non hispanic white
whiteids<-data_gail_final_white$id
data_gail_final_black<-data_gail_final2[which (data_gail_final2$Race==2),] #CARE model
blackids<-data_gail_final_black$id
data_gail_final_asian<-data_gail_final2[which (data_gail_final2$Race==6),] #AA model
asianids<-data_gail_final_asian$id


#Plot distributions
# #Plot relative risk distribution and absolute risk distribution by race
# #non hispanic white
# ggplot(data_gail_final_white, aes(x=Relrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Relative Risk",y="Frequency") + ggtitle("Relative Risk Distribution for Gail Model Non-Hispanic White")
# ggsave("filepath/Gail_relrisk_white.png", width=4, height=4)
# 
# ggplot(data_gail_final_white, aes(x=absrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Absolute Risk",y="Frequency") + ggtitle("Absolute Risk Distribution for Gail Model Non-Hispanic White")
# ggsave("filepath/Gail_absrisk_white.png", width=4, height=4)
# 
# #Among black americans
# ggplot(data_gail_final_black, aes(x=Relrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Relative Risk",y="Frequency") + ggtitle("Relative Risk Distribution for CARE Model")
# ggsave("filepath/Gail_relrisk_black.png", width=4, height=4)
# 
# ggplot(data_gail_final_black, aes(x=absrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Absolute Risk",y="Frequency") + ggtitle("Absolute Risk Distribution for CARE Model")
# ggsave("filepath/Gail_absrisk_black.png", width=4, height=4)


# #Among asian americans
# ggplot(data_gail_final_asian, aes(x=Relrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Relative Risk",y="Frequency") + ggtitle("Relative Risk Distribution for Asian American Gail Model")
# ggsave("filepath/Gail_relrisk_asian.png", width=4, height=4)
# 
# ggplot(data_gail_final_asian, aes(x=absrisk_gail))+ geom_histogram(binwidth=0.5, fill="lightblue")+labs(x="Absolute Risk",y="Frequency") + ggtitle("Absolute Risk Distribution for Asian American Gail Model")
# ggsave("filepath/Gail_absrisk_asian.png", width=4, height=4)


########################################################
# Calibration of relative risk - use iCARE Validate    #
########################################################

#Use data2 -assign values to align with iCARE 
#assume all follow up covers 5 years

#Need time of onset added to dataset

data_validategail<-data2[,c("id","cohort","case","race","timediff_case_yr")]

#If timediff case yr is missing assign as infinite (past 5 y)

data_validategail$time.of.onset<-ifelse(data_validategail$case==1, data_validategail$timediff_case_yr, Inf)

    #Data check!
    #table(data_validategail$case, data_validategail$time.of.onset,useNA="always")
    #table(data2$timediff_case_yr, data2$case_5y, useNA="always")

#Combine with relative and absolute risk data
#head(data_validategail)
#For validation cannot have NAs - delete if time.of.onset==NA

data_validategail_full<-data_validategail[which (!is.na(data_validategail$time.of.onset)),]
delete.ids<-setdiff(data_validategail$id, data_validategail_full$id)

data_gail_final2new<-data_gail_final2[which (! data_gail_final2$id %in% delete.ids),]

data_validategail_white<-data_validategail_full[which (data_validategail_full$id %in% whiteids & !data_validategail_full$id %in% delete.ids),]
data_validategail_black<-data_validategail_full[which (data_validategail_full$id %in% blackids & !data_validategail_full$id %in% delete.ids),]
data_validategail_asian<-data_validategail_full[which (data_validategail_full$id %in% asianids & !data_validategail_full$id %in% delete.ids),]

data_gail_final_whitenew<-data_gail_final2new[which (data_gail_final2new$id %in% whiteids),]
data_gail_final_blacknew<-data_gail_final2new[which (data_gail_final2new$id %in% blackids),]
data_gail_final_asiannew<-data_gail_final2new[which (data_gail_final2new$id %in% asianids),]

#########
#Samples

data_validategail_full<-cbind(data_validategail_full, data_gail_final2new)
data_validategail_full$sampling.weights<-1

#Limit to races

data_validategail_white<-cbind(data_validategail_white, data_gail_final_whitenew)
data_validategail_white$sampling.weights<-1

data_validategail_black<-cbind(data_validategail_black, data_gail_final_blacknew)
data_validategail_black$sampling.weights<-1

data_validategail_asian<-cbind(data_validategail_asian, data_gail_final_asiannew)
data_validategail_asian$sampling.weights<-1

#Rename variables 
data_validategail_full<-data_validategail_full[,c("id","T1","T2","case","Relrisk_gail","absrisk_gailnew","time.of.onset","sampling.weights")]
data_validategail_white<-data_validategail_white[,c("id","T1","T2","case","Relrisk_gail","absrisk_gailnew","time.of.onset","sampling.weights")]
data_validategail_black<-data_validategail_black[,c("id","T1","T2","case","Relrisk_gail","absrisk_gailnew","time.of.onset","sampling.weights")]
data_validategail_asian<-data_validategail_asian[,c("id","T1","T2","case","Relrisk_gail","absrisk_gailnew","time.of.onset","sampling.weights")]

#Rename columns to fit iCARE scheme:
#abs risk given as % in gail data

colnames(data_validategail_full)<-c("id","study.entry.age","study.exit.age","observed.outcome","Relrisk_gail","absrisk_gailnew", "time.of.onset","sampling.weights")
#Subset to needed vars
data_validategail_full2<-data_validategail_full[,c("study.entry.age","study.exit.age","observed.outcome","time.of.onset","sampling.weights")]

    #Data check!
    # mean(data_validategail_full$observed.outcome, na.rm=TRUE) #overall
    # mean(data_validategail_full$absrisk_gailnew, na.rm=TRUE) #5 year


colnames(data_validategail_white)<-c("id","study.entry.age","study.exit.age","observed.outcome","Relrisk_gail","absrisk_gailnew", "time.of.onset","sampling.weights")
data_validategail_white2<-data_validategail_white[,c("study.entry.age","study.exit.age","observed.outcome","time.of.onset","sampling.weights")]

colnames(data_validategail_black)<-c("id","study.entry.age","study.exit.age","observed.outcome","Relrisk_gail","absrisk_gailnew", "time.of.onset","sampling.weights")
data_validategail_black2<-data_validategail_black[,c("study.entry.age","study.exit.age","observed.outcome","time.of.onset","sampling.weights")]

colnames(data_validategail_asian)<-c("id","study.entry.age","study.exit.age","observed.outcome","Relrisk_gail","absrisk_gailnew", "time.of.onset","sampling.weights")
data_validategail_asian2<-data_validategail_asian[,c("study.entry.age","study.exit.age","observed.outcome","time.of.onset","sampling.weights")]

#Vector risk assignments
#Linear predict assignments

gailrisk=as.vector(data_validategail_full$absrisk_gailnew)
linearpredict=as.vector(log(data_validategail_full$Relrisk_gail))
gailriskwhite=as.vector(data_validategail_white$absrisk_gailnew)
linearpredictwhite=as.vector(log(data_validategail_white$Relrisk_gail))
gailriskblack=as.vector(data_validategail_black$absrisk_gailnew)
linearpredictblack=as.vector(log(data_validategail_black$Relrisk_gail))
gailriskasian=as.vector(data_validategail_asian$absrisk_gailnew)
linearpredictasian=as.vector(log(data_validategail_asian$Relrisk_gail))


#####################################################
# Apply iCARE validate functions to calibrate model
# This is based on deciles of relative risk, lps
#####################################################


#Will only run if sampling.weights set to NULL

data_validategail_full2$sampling.weights<-NULL

#Result for ALL races, full cohort
Gail_AllRace_Result=ModelValidation(data_validategail_full2, 
                                    total.followup.validation = FALSE,
                                    predicted.risk = gailrisk, #Vector of absolute risk, supplied
                                    predicted.risk.interval = 5, 
                                    linear.predictor = linearpredict, 
                                    iCARE.model.object = 
                                      list(model.formula = NULL,
                                           model.cov.info = NULL,
                                           model.snp.info = NULL,
                                           model.log.RR = NULL,
                                           model.ref.dataset = NULL,
                                           model.ref.dataset.weights = NULL,
                                           model.disease.incidence.rates = NULL,
                                           model.competing.incidence.rates = NULL,
                                           model.bin.fh.name = NA,
                                           apply.cov.profile  = NULL,
                                           apply.snp.profile = NULL, 
                                           n.imp = NULL, use.c.code = NULL,
                                           return.lp = TRUE, 
                                           return.refs.risk = TRUE),
                                    number.of.percentiles = 10,
                                    reference.entry.age = NULL, 
                                    reference.exit.age = NULL,
                                    predicted.risk.ref = NULL,
                                    linear.predictor.ref = NULL,
                                    linear.predictor.cutoffs = NULL,
                                    dataset = "Gail All Dataset", 
                                    model.name = "Gail All Prediction Model")   


#Prints Summary Calibration  
Gail_AllRace_Result
capture.output(Gail_AllRace_Result, file="filepath/Gail_AllRace_Calibration.txt",sep=',')

#Gives expected and observed by decile
Gail_AllRace_Result$Category_Results



#Plot the calibration results
par(mar=c(1,1,1,1))
#Note: Replace "STUDY" with the acronym for your study
#This uses provided plotModelValidation code with iCARE, in race-specific versions use modified function plotModelValidation_withNAs provided in this Github folder

Gail_AllRace_CalibrationPlot_STUDY=plotModelValidation(data_validategail_full2, Gail_AllRace_Result,
                                                         dataset="STUDY",
                                                         model.name="Gail all races",
                                                         x.lim.absrisk = NULL,
                                                         y.lim.absrisk = NULL, 
                                                         x.lab.absrisk = "Expected Absolute Risk (%)", 
                                                         y.lab.absrisk = "Observed Absolute Risk (%)", 
                                                         x.lim.RR = NULL,
                                                         y.lim.RR = NULL, 
                                                         x.lab.RR = "Expected Relative Risk", 
                                                         y.lab.RR = "Observed Relative Risk",
                                                         risk.score.plot.kernel = "gaussian",
                                                         risk.score.plot.bandwidth = "nrd0",
                                                         risk.score.plot.percent.smooth = 50)



############################################################################################################
# Repeat for White alone, Black alone (CARE model), and Asian American alone (Asian-American Risk Model)
#############################################################################################################

#White alone

data_validategail_white2$sampling.weights<-NULL

Gail_WhiteRace_Result=ModelValidation(data_validategail_white2, #Name of dataset changes - white only
                                      total.followup.validation = FALSE,
                                      predicted.risk = gailriskwhite, #Vector of absolute risk, supplied, changes for white only
                                      predicted.risk.interval = 5, 
                                      linear.predictor = linearpredictwhite, #This changes specific to input dataset
                                      iCARE.model.object = 
                                        list(model.formula = NULL,
                                             model.cov.info = NULL,
                                             model.snp.info = NULL,
                                             model.log.RR = NULL,
                                             model.ref.dataset = NULL,
                                             model.ref.dataset.weights = NULL,
                                             model.disease.incidence.rates = NULL,
                                             model.competing.incidence.rates = NULL,
                                             model.bin.fh.name = NA,
                                             apply.cov.profile  = NULL,
                                             apply.snp.profile = NULL, 
                                             n.imp = NULL, use.c.code = NULL,
                                             return.lp = TRUE, 
                                             return.refs.risk = TRUE),
                                      number.of.percentiles = 10,
                                      reference.entry.age = NULL, 
                                      reference.exit.age = NULL,
                                      predicted.risk.ref = NULL,
                                      linear.predictor.ref = NULL,
                                      linear.predictor.cutoffs = NULL,
                                      dataset = "Gail All Dataset", 
                                      model.name = "Gail All Prediction Model")   


#IF ALL PLOTS do not show then NaNs are likely in results 
#Function to plot by replacing Nas with 0- use modified function "plotModelValidationwithNAs"

#Download function from github: 
source("filepath/PlotCalibration_withNAs.R")

#All NaNs must be changed to "0" for plotting to occur
#Create new "Category results"
Category.Results<-Gail_WhiteRace_Result$Category_Results
Category.Results[Category.Results=="NaN"]<-0

Gail_WhiteRace_CalibrationPlot_STUDY=plotModelValidationwithNAs(data_validategail_white2, Category.Results, Gail_WhiteRace_Result,
                                                                  dataset="STUDY", #Name as cohort
                                                                  model.name="Gail - White",
                                                                  x.lim.absrisk = NULL,
                                                                  y.lim.absrisk = NULL, 
                                                                  x.lab.absrisk = "Expected Absolute Risk (%)", 
                                                                  y.lab.absrisk = "Observed Absolute Risk (%)", 
                                                                  x.lim.RR = NULL,
                                                                  y.lim.RR = NULL, 
                                                                  x.lab.RR = "Expected Relative Risk", 
                                                                  y.lab.RR = "Observed Relative Risk",
                                                                  risk.score.plot.kernel = "gaussian",
                                                                  risk.score.plot.bandwidth = "nrd0",
                                                                  risk.score.plot.percent.smooth = 50)

Gail_WhiteRace_Result
capture.output(Gail_WhiteRace_Result, file="filepath/Gail_WhiteRace_Calibration.txt",sep=',')

#Black alone

data_validategail_black2$sampling.weights<-NULL

Gail_blackRace_Result=ModelValidation(data_validategail_black2, #Name of dataset changes - black only
                                      total.followup.validation = FALSE,
                                      predicted.risk = gailriskblack, #Vector of absolute risk, supplied, changes for black only
                                      predicted.risk.interval = 5, 
                                      linear.predictor = linearpredictblack, #This changes specific to input dataset
                                      iCARE.model.object = 
                                        list(model.formula = NULL,
                                             model.cov.info = NULL,
                                             model.snp.info = NULL,
                                             model.log.RR = NULL,
                                             model.ref.dataset = NULL,
                                             model.ref.dataset.weights = NULL,
                                             model.disease.incidence.rates = NULL,
                                             model.competing.incidence.rates = NULL,
                                             model.bin.fh.name = NA,
                                             apply.cov.profile  = NULL,
                                             apply.snp.profile = NULL, 
                                             n.imp = NULL, use.c.code = NULL,
                                             return.lp = TRUE, 
                                             return.refs.risk = TRUE),
                                      number.of.percentiles = 10,
                                      reference.entry.age = NULL, 
                                      reference.exit.age = NULL,
                                      predicted.risk.ref = NULL,
                                      linear.predictor.ref = NULL,
                                      linear.predictor.cutoffs = NULL,
                                      dataset = "Gail All Dataset", 
                                      model.name = "Gail All Prediction Model")   


#Prints Summary Calibration  
Gail_blackRace_Result
capture.output(Gail_blackRace_Result, file="filepath/Gail_blackRace_Calibration.txt",sep=',')

#Gives expected and observed by decile
#Gail_blackRace_Result$Category_Results

#Plot the calibration results
par(mar=c(1,1,1,1))

#Note: Replace "STUDY" with the acronym for your study

Category.Results<-Gail_blackRace_Result$Category_Results
Category.Results[Category.Results=="NaN"]<-0

Gail_blackRace_CalibrationPlot_STUDY=plotModelValidationwithNAs(data_validategail_black2, Category.Results,Gail_blackRace_Result,
                                                                  dataset="STUDY", #Name as cohort
                                                                  model.name="Gail - black",
                                                                  x.lim.absrisk = NULL,
                                                                  y.lim.absrisk = NULL, 
                                                                  x.lab.absrisk = "Expected Absolute Risk (%)", 
                                                                  y.lab.absrisk = "Observed Absolute Risk (%)", 
                                                                  x.lim.RR = NULL,
                                                                  y.lim.RR = NULL, 
                                                                  x.lab.RR = "Expected Relative Risk", 
                                                                  y.lab.RR = "Observed Relative Risk",
                                                                  risk.score.plot.kernel = "gaussian",
                                                                  risk.score.plot.bandwidth = "nrd0",
                                                                  risk.score.plot.percent.smooth = 50)



#Asian alone

data_validategail_asian2$sampling.weights<-NULL

Gail_asianRace_Result=ModelValidation(data_validategail_asian2, #Name of dataset changes - asian only
                                      total.followup.validation = FALSE,
                                      predicted.risk = gailriskasian, #Vector of absolute risk, supplied, changes for asian only
                                      predicted.risk.interval = 5, 
                                      linear.predictor = linearpredictasian, #This changes specific to input dataset
                                      iCARE.model.object = 
                                        list(model.formula = NULL,
                                             model.cov.info = NULL,
                                             model.snp.info = NULL,
                                             model.log.RR = NULL,
                                             model.ref.dataset = NULL,
                                             model.ref.dataset.weights = NULL,
                                             model.disease.incidence.rates = NULL,
                                             model.competing.incidence.rates = NULL,
                                             model.bin.fh.name = NA,
                                             apply.cov.profile  = NULL,
                                             apply.snp.profile = NULL, 
                                             n.imp = NULL, use.c.code = NULL,
                                             return.lp = TRUE, 
                                             return.refs.risk = TRUE),
                                      number.of.percentiles = 10,
                                      reference.entry.age = NULL, 
                                      reference.exit.age = NULL,
                                      predicted.risk.ref = NULL,
                                      linear.predictor.ref = NULL,
                                      linear.predictor.cutoffs = NULL,
                                      dataset = "Gail All Dataset", 
                                      model.name = "Gail All Prediction Model")   


#Prints Summary Calibration  
Gail_asianRace_Result
capture.output(Gail_asianRace_Result, file="filepath/Gail_AsianRace_Calibration.txt",sep=',')

#Gives expected and observed by decile
Gail_asianRace_Result$Category_Results

#Plot the calibration results
par(mar=c(1,1,1,1))
#Note: Replace "STUDY" with the acronym for your study
Category.Results<-Gail_asianRace_Result$Category_Results
Category.Results[Category.Results=="NaN"]<-0

Gail_asianRace_CalibrationPlot_STUDY=plotModelValidation(data_validategail_asian2, Category.Results,Gail_asianRace_Result,
                                                           dataset="STUDY", #Name as cohort
                                                           model.name="Gail - asian",
                                                           x.lim.absrisk = NULL,
                                                           y.lim.absrisk = NULL, 
                                                           x.lab.absrisk = "Expected Absolute Risk (%)", 
                                                           y.lab.absrisk = "Observed Absolute Risk (%)", 
                                                           x.lim.RR = NULL,
                                                           y.lim.RR = NULL, 
                                                           x.lab.RR = "Expected Relative Risk", 
                                                           y.lab.RR = "Observed Relative Risk",
                                                           risk.score.plot.kernel = "gaussian",
                                                           risk.score.plot.bandwidth = "nrd0",
                                                           risk.score.plot.percent.smooth = 50)

