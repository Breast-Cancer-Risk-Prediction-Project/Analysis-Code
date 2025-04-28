#External code - additional data cleaning to allow for validation 
#This code will create the validation dataset with appropriate variables and formats
#and participant exclusions (e.g., >=20y, <=75y, only women) that can then be used in
#the calibration codes for each model.

#  Notes: Please fill in "FILEPATH" where indicated and in Step 1 please enter the cohort acronym for
# "ACRONYM"


rm(list=ls())
#Required packages/libraries
library('tidyverse')
library('dplyr')
library("lubridate")

###############################################################################  
#1. Read in your dataset - data_core and data_brca
###############################################################################  

data_core<-read.csv("/FILEPATH.csv") 
data_brca<-read.csv("/FILEPATH.csv")

#PLEASE INSERT YOUR COHORT ACRONYM HERE

data_core$cohort<-"ACRONYM"

###############################################################################  
#2. You must ensure ID formats match for core and brca datasets. 
############################################################################### 

#Format for subject_id=ACRONYM_ID, e.g.: "NHS_X0021190"
#Format for id= ID,                e.g.: "X0021190"

#Ensure class is character
if (class(data_core$subject_id) !="character"){print ("Error: Subject ID not character")}

  #Note: If this error occurs please check your data

#Data checks - ensure IDs align
  
  if("subject_id" %in% colnames(data_brca)){ matchvar="subject_id" }
  if("id" %in% colnames(data_brca)){ matchvar="id" }
  if("subject_id" %in% colnames(data_brca) & "id" %in% colnames(data_brca)){matchvar="both"}
  
  #Change id to character if numeric
  if (class(data_core$id) !="character"){data_core$subject_id<-as.character(data_core$id)}

  #The following code will warn you if the formats do not match between CORE and BRCA

  if (matchvar=="id"|matchvar=="both"){
    if (nchar(data_core$id[1]) != nchar(data_brca$id[1])){print("Error, core and brca IDs do not align")}
  }
  if (matchvar=="subject_id"|matchvar=="both"){
    if (nchar(data_core$subject_id[1]) != nchar(data_brca$subject_id[1])){print("Error, core and brca subject IDs do not align")}
  }

#If these errors print please fix formatting of ID variables before moving on. 

#Coding solutions to common mismatches are included below. 
#The datasets to introduce changes to will depend on any errors.

  #1.Lack of leading zeros in the ID variable (e.g. subject_id="COHORT_000213 but id=213)
  
  #Here "N" is the number of leading zeros that need to be added - e.g. "0006", N=4 -- change accordingly
  
  #  data$newid<-str_pad(data$id, N, side="left", pad="0")
  #  drop<-c("id")
  #  data<-data[,!(names(data) %in% drop),]
  #  data<-data %>%
  #    rename(id=newid)

  #2. Missing "COHORT_" leader in subject_id
  #  data$subject_id<-paste0("COHORT","_", data$id)

  #3. Included "-" instead of "_":
  # data$subject_id<-gsub("-", "_", data$subject_id)

############################################################################### 

###############################################################################    
#3.Combine core and incident brca datasets
#  This will choose the correct merge function for your data structure. 
#  Option 1 is chosen if you have ALL ids in the BRCA dataset. Option 2 is 
#  Chosen if you have only CASES in the BRCA dataset.
###############################################################################    

  #Merge functions:
  casemerge1a<-function(data1,data2, varlist){
    data1$dxdatemiss<-ifelse(data1$dxdate_primary1 ==8888|data1$dxdate_primary1 ==888 | data1$dxdate_primary1 ==7777| data1$dxdate_primary1 ==777| is.na(data1$dxdate_primary1), 1,0)
    data1$stagemiss<-ifelse(is.na(data1$stage_primary1) | data1$stage_primary1==777 | data1$stage_primary1==888, 1, 0)
    data1$detectionmiss<-ifelse(is.na(data1$detection_primary1)| data1$detection_primary1==777 | data1$detection_primary1==888, 1, 0)
    cases<-data1[which ((data1$dxdatemiss==0) | (data1$dxdatemiss==1 & data1$stagemiss==0 | data1$detectionmiss==0)),]
    caselist<-cases$id
    data1$case<-ifelse(data1$id %in% caselist, 1, 0)
    case1<-data1[which (data1$id %in% caselist),]
    case2<-data2[which (data2$id %in% caselist),]
    df<-merge(case2, case1, by=varlist)
    names(df)<-tolower(names(df))
    return(df)}
  casemerge1b<-function(data1,data2, varlist){
    data1$dxdatemiss<-ifelse(data1$dxdate_primary1 ==8888|data1$dxdate_primary1 ==888 | data1$dxdate_primary1 ==7777| data1$dxdate_primary1 ==777| is.na(data1$dxdate_primary1), 1,0)
    data1$stagemiss<-ifelse(is.na(data1$stage_primary1) | data1$stage_primary1==777 | data1$stage_primary1==888, 1, 0)
    data1$detectionmiss<-ifelse(is.na(data1$detection_primary1)| data1$detection_primary1==777 | data1$detection_primary1==888, 1, 0)
    cases<-data1[which ((data1$dxdatemiss==0) | (data1$dxdatemiss==1 & data1$stagemiss==0 | data1$detectionmiss==0)),]
    caselist<-cases$subject_id
    data1$case<-ifelse(data1$subject_id %in% caselist, 1, 0)
    case1<-data1[which (data1$subject_id %in% caselist),]
    case2<-data2[which (data2$subject_id %in% caselist),]
    df<-merge(case2, case1, by=varlist)
    names(df)<-tolower(names(df))
    return(df)}
  casemerge2a<-function(data1,data2, varlist){
    cases<-data1
    caselist<-cases$id
    data1$case<-ifelse(data1$id %in% caselist, 1, 0)
    case1<-data1[which (data1$id %in% caselist),]
    case2<-data2[which (data2$id %in% caselist),]
    df<-merge(case2, case1, by=varlist)
    names(df)<-tolower(names(df))
    return(df)}
  casemerge2b<-function(data1,data2, varlist){
    cases<-data1
    caselist<-cases$subject_id
    data1$case<-ifelse(data1$subject_id %in% caselist, 1, 0)
    case1<-data1[which (data1$subject_id %in% caselist),]
    case2<-data2[which (data2$subject_id %in% caselist),]
    df<-merge(case2, case1, by=varlist)
    names(df)<-tolower(names(df))
    return(df)}

  controlmerge1a<-function(data1,data2,data3,varlist){
    caselist<-data2$id
    controls<-data1[which (!(data1$id) %in% caselist),]
    controllist<-controls$id
    control1<-data1[which (data1$id %in% controllist),]
    control2<-data3[which (data3$id %in% controllist),]
    df<-merge(control1, control2, by=varlist)
    df$case<-0
    df$dxdatemiss<-NA
    df$stagemiss<-NA
    df$detectionmiss<-NA
    names(df)<-tolower(names(df))
    return(df)
  }
  controlmerge1b<-function(data1,data2,data3,varlist){
    caselist<-data2$subject_id
    controls<-data1[which (!(data1$subject_id) %in% caselist),]
    controllist<-controls$subject_id
    control1<-data1[which (data1$subject_id %in% controllist),]
    control2<-data3[which (data3$subject_id %in% controllist),]
    df<-merge(control1, control2, by=varlist)
    df$case<-0
    df$dxdatemiss<-NA
    df$stagemiss<-NA
    df$detectionmiss<-NA
    names(df)<-tolower(names(df))
    return(df)
  }
  controlmerge2a<-function(data1,data2, varlist){
    a<-setdiff(data1$id, data2$id)
    data_allcontrol<-data1[which (data1$id %in% (a)),]
    controllist<-data_allcontrol$id
    df<-data1[which (data1$id %in% controllist),]
    df$case<-0
    names(df)<-tolower(names(df))
    return(df)
  }
  controlmerge2b<-function(data1,data2,varlist){
    a<-setdiff(data1$subject_id, data2$subject_id)
    data_allcontrol<-data1[which (data1$subject_id %in% (a)),]
    controllist<-data_allcontrol$subject_id
    df<-data1[which (data1$subject_id %in% controllist),]
    df$case<-0
    names(df)<-tolower(names(df))
    return(df)
  }

#3a. Merge controls and cases together step 1
if (dim(data_brca)[1]==dim(data_core)[1]){
  if (matchvar=="both"){   #cases and controls are both in BRCA
          data_allcase<-casemerge1a(data_brca, data_core, c("subject_id","id"))}
  else if (matchvar=="id"){
    data_allcase<-casemerge1a(data_brca, data_core, c("id"))
  }
  else if (matchvar=="subject_id"){
    data_allcase<-casemerge1b(data_brca, data_core, c("subject_id"))
  }
}      

  if (dim(data_brca)[1]==dim(data_core)[1]){
    if (matchvar=="both"){   #cases and controls are both in BRCA
      data_allcontrol<-controlmerge1a(data_core,data_allcase,data_brca, c("subject_id","id"))}
    else if (matchvar=="id"){
      data_allcontrol<-controlmerge1a(data_core,data_allcase,data_brca, c("id"))
    }
    else if (matchvar=="subject_id"){
      data_allcontrol<-controlmerge1b(data_core,data_allcase,data_brca, c("subject_id"))
    }
  }  

#BRCA only had cases to start
  
  if (dim(data_brca)[1]<dim(data_core)[1]){
    if (matchvar=="both"){   #cases and controls are both in BRCA
      data_allcase<-casemerge2a(data_brca,data_core, c("subject_id","id"))}
    else if (matchvar=="id"){
      data_allcase<-casemerge2a(data_brca,data_core, c("id"))
    }
    else if (matchvar=="subject_id"){
      data_allcase<-casemerge2b(data_brca, data_core,c("subject_id"))
    }
  }      

  if (dim(data_brca)[1]<dim(data_core)[1]){
    if (matchvar=="both"){   #cases and controls are both in BRCA
      data_allcontrol<-controlmerge2a(data_core,data_brca, c("subject_id","id"))}
    else if (matchvar=="id"){
      data_allcontrol<-controlmerge2a(data_core,data_brca, c("id"))
    }
    else if (matchvar=="subject_id"){
      data_allcontrol<-controlmerge2b(data_core,data_brca, c("subject_id"))
    }
  }  

  
#3b. Merge controls and cases -step 2
# Make sure there are no differences in case and control variables. If so, will need to add before appending

#Note: this will add case variables as NA for controls, if ONLY cases were in original brca dataset
  if(dim(data_brca)[1]<dim(data_core)[1]){
    a<-setdiff(names(data_brca),names(data_core))
    for (var in a){
    data_allcontrol[,var]<-NA
  }
  }

  
  #Set as DF
  data_allcase<-as.data.frame(data_allcase)
  data_allcontrol<-as.data.frame(data_allcontrol)
  

  #3c. Make sure ALL variable names in DD are included in dataset- all lower case from earlier commands
  names1<-c("subject_id","id","record_date","qcycle","baseline","lastfup","lastfup_reason",
            "status","birth_year","sex","age","race","ethnicity","education","ajancestry",
            "height","weight","bmi","waist","hip","whr","bmi_earlyadult","alcohol_status",
            "alcohol_init","alcohol_amt","alcohol_stop","alcohol_dur","smoking_status",
            "smoking_init","smoking_amt","smoking_dur","smoking_stop","fhx_fdr_brca","sisters",
            "brcancersis","daughters","brcancerdau","brcancermom","brcancerdad","fhisfstbc",
            "fhisfstbcnr","fam1grbc50","ovcancersis","ovcancerdau","ovcancermom","fhisfstoc",
            "fhisfstocnr","fam1groc50","biopsies_yesno","biopsies_number","bbd_history",
            "bbd_number","bbd_type1","bbd_year1","bbd_type2","bbd_year2","bbd_type3","bbd_year3",
            "bbd_type4","bbd_year4","agemenarche","parous","age_preg1","parity","age_preg2",
            "age_preg3","age_preg4","age_preg5","age_preg6","age_preg7","age_preg8",
            "age_preg9","age_preg10","breastfeed","breastfeed_dur","breastfeed_dur_b1",
            "breastfeed_dur_b2","breastfeed_dur_b3","breastfeed_dur_b4","breastfeed_dur_b5",
            "breastfeed_dur_b6","breastfeed_dur_b7","breastfeed_dur_b8","breastfeed_dur_b9",
            "breastfeed_dur_b10","ocuse_ever","ocuse_current","ocuse_dur","ocuse_start","ocuse_stop",
            "othcontracep_ever","othcontracep_current","meno_status","meno_age","meno_reason",
            "hrtuse","hrt_dur","hrtuse_ep","hrtep_dur","hrtuse_eonly","hrteonly_dur","pa_mets","pa_pct",
            "screen_ever","screen_start","lastscreen_year","dxdate_primary1","detection_primary1",
            "detection_detail_primary1","invasive_primary1","dxdate_primary2","detection_primary2",
            "detection_detail_primary2","invasive_primary2","laterality_primary2","lateralitytime_primary2",
            "stage_primary1","grade_primary1","size_primary1","sizecat_primary1","er_primary1",
            "pr_primary1","her2_primary1","ki67_primary1","ki67cat_primary1","stage_primary2",
            "grade_primary2","size_primary2","sizecat_primary2","er_primary2","pr_primary2",
            "her2_primary2","ki67_primary2","ki67cat_primary2")
  
  
  #Add columns if they aren't in dataset
  
  diff<-setdiff(names1, colnames(data_allcase))
  for (var in diff){
    data_allcase[,var]<-NA
  }
  
  diff
  
  diff<-setdiff(names1, colnames(data_allcontrol))
  for (var in diff){
    data_allcontrol[,var]<-NA
  }
  
  #3d. Final merge
  data<-rbind(data_allcontrol, data_allcase)
  
  #End of step 1 data cleaning
        
        #Step 2 data cleaning:
#######################################################
#4. Date format alignment - missing to NAs - and calculations
#   Note: all date formats must match DD for this code to work.
#######################################################
    
  data<-data %>%
    
    #record_year will only be used here
    
    mutate(record_date=ifelse(record_date=="08/08/8000", NA, record_date),
           record_year=substr(record_date, 7,10), 
           record_year=as.numeric(record_year),
           
           #Change dx year so it is year not days for WHI
           dxyear_primary1=ifelse(cohort=="WHI" & (dxdate_primary1==8888| dxdate_primary1==888 | dxdate_primary1==7777 | dxdate_primary1==777 ), NA, round((1900+dxdate_primary1/365),0)),
           dxyear_primary2=ifelse(cohort=="WHI" & (dxdate_primary2==8888| dxdate_primary2==888 | dxdate_primary2==7777 | dxdate_primary2==777), NA, round((1900+dxdate_primary2/365),0)),
           
           #change variable to dxyear, as numeric for all other cohorts
           dxyear_primary1=ifelse(dxdate_primary1==8888| dxdate_primary1==888 | dxdate_primary1==7777 | dxdate_primary1==777, NA, suppressWarnings(as.numeric(dxdate_primary1))),
           
           #label as prevalent case if dx is <= record_year
           prevalent_case_primary1=ifelse(!is.na(dxyear_primary1) & !is.na(record_year) & dxyear_primary1<=record_year,1, 0),
           
           #Repeat for 2nd dx
           dxyear_primary2=ifelse(dxdate_primary2==8888| dxdate_primary2==888 | dxdate_primary2==7777 | dxdate_primary2==777, NA, suppressWarnings(as.numeric(dxdate_primary2))),
           prevalent_case_primary2=ifelse(!is.na(dxyear_primary2) & !is.na(record_year) & dxyear_primary2<=record_year,1, 0),
           #change birth_year, bbd year, lastscreen year, lastfup year to NA if marked as missing
           birth_year=ifelse(birth_year==8888, NA, as.numeric(birth_year)),
           bbd_year1=ifelse(bbd_year1==8888| bbd_year1==7777| bbd_year1==888, NA, bbd_year1),
           bbd_year2=ifelse(bbd_year2==8888| bbd_year2==7777| bbd_year2==888, NA, bbd_year2),
           bbd_year3=ifelse(bbd_year3==8888| bbd_year3==7777| bbd_year3==888, NA, bbd_year3),
           bbd_year4=ifelse(bbd_year4==8888| bbd_year4==7777| bbd_year4==888, NA, bbd_year4),
           lastscreen_year=ifelse(lastscreen_year==8888| lastscreen_year==7777| lastscreen_year==888| lastscreen_year==777, NA, lastscreen_year),
           lastfup=ifelse(lastfup>=7777, NA, lastfup),
           
           #Calculate time between enroll and dx by year only - 
           #Will use this variable to exclude individuals who were diagnosed within 1 year
           time_bltodx1=dxyear_primary1-record_year,
           time_bltodx2=dxyear_primary2-record_year,
           #Fix age if missing but record year and birth year are available
           age=ifelse(age==888, NA, age),
           age=ifelse(is.na(age) & !is.na(record_year) & !is.na(birth_year), record_year-birth_year, age ))
  ########################################################
  #5. Ensure certain variables are numeric, if not already
  ########################################################
  
  
  data<-as.data.frame(data)
  
  asnumeric<-function(data, varlist){
    for (i in varlist){
      data[,i]<-as.numeric(data[,i])
    }
    return(data)
  }
  
  varlist1<-c("height","weight","bmi","waist","hip","whr","bmi_earlyadult","alcohol_init","alcohol_amt","alcohol_stop",
              "alcohol_dur","smoking_dur","smoking_amt","smoking_stop","smoking_init","agemenarche","parity","parous",
              "age_preg1","age_preg2","age_preg3","age_preg4","age_preg5","age_preg6","age_preg7","age_preg8","age_preg9","age_preg10",
              "breastfeed_dur_b1","breastfeed_dur_b2","breastfeed_dur_b3","breastfeed_dur_b4","breastfeed_dur_b5","breastfeed_dur_b6","breastfeed_dur_b7","breastfeed_dur_b8","breastfeed_dur_b9","breastfeed_dur_b10",
              "ocuse_dur","ocuse_start","ocuse_stop","hrt_dur","hrtep_dur","hrteonly_dur","pa_mets","pa_pct",
              "size_primary1","size_primary2","ki67_primary1","ki67_primary2")
  
  data<-suppressWarnings(asnumeric(data, varlist1))
  
  ###################################################
  #6. Outlier assignment 
  ################################################### 
  
  data1<- data %>%
    mutate(
      age_outlier    =ifelse(age<18 & !is.na(age), 1,
                             ifelse(!is.na(age) & age>100, 1, 0)),
      height_outlier        =ifelse(height<120, 1,
                                    ifelse(height<666 & height>200, 1, 0)),
      weight_outlier        =ifelse(weight<35,1,
                                    ifelse(weight<666 & weight>150, 1, 0)),
      bmi_outlier           =ifelse(bmi<=15, 1,
                                    ifelse(bmi<666& bmi>50, 1, 0)),
      waist_outlier         =ifelse(waist<=50, 1,
                                    ifelse(waist<666& waist>160, 1, 0)),
      hip_outlier           =ifelse(hip<70, 1,
                                    ifelse(hip<666 & hip>160, 1, 0)),
      wtohip_outlier        =ifelse(whr<=0,1,
                                    ifelse(whr<666 & whr>1.2, 1, 0)),
      bmiearly_outlier      =ifelse(bmi_earlyadult<=0, 1, 
                                    ifelse(bmi_earlyadult>50 & bmi_earlyadult<666, 1, 0)),
      #changed initiation to >90, <10 as outliers
      alcoholinit_outlier   =ifelse(alcohol_init<10, 1,
                                    ifelse(alcohol_init<666 & alcohol_init>90, 1, 0)),
      alcoholamt_outlier    =ifelse(alcohol_amt<0, 1,
                                    ifelse(alcohol_amt<666 & alcohol_amt>50, 1, 0)),
      #changed this to >50 since the 7 category variable allows up to 45+
      #changed alcohol stop to <10
      alcoholstop_outlier   =ifelse(alcohol_stop<10,1,
                                    ifelse(alcohol_stop<666 & alcohol_stop>100, 1, 0)),
      #changed duration outlier to fit age
      alcoholdur_outlier    =ifelse(alcohol_dur<0, 1, 
                                    ifelse(alcohol_dur<666 & alcohol_dur>(age-10), 1, 0)),
      #added outlier here for CPS3 - assuming 999 are missing
      #smoking dur changed to reflect age and assuming smoking can start at age 10 as opposed to age 18
      smokingdur_outlier    =ifelse(smoking_dur<0, 1, 
                                    ifelse(smoking_dur<666 & smoking_dur>(age-10), 1,
                                           ifelse(smoking_dur>999, 1, 0))),
      smokingamt_outlier    =ifelse(smoking_amt<0, 1,
                                    ifelse(smoking_amt<666 & smoking_amt>60, 1, 0)),
      #change smoking stop
      smokingstop_outlier   =ifelse(smoking_stop<10, 1,
                                    ifelse(smoking_stop<666 & smoking_stop>100, 1, 0)),
      #change to 10 
      smokinginit_outlier   =ifelse(smoking_init<10, 1,
                                    ifelse(smoking_init<666 & smoking_init>90, 1, 0)),
      
      agemenarche_outlier   =ifelse(agemenarche<8, 1,
                                    ifelse(agemenarche<666 & agemenarche>18, 1, 0)),
      parity_outlier        =ifelse(parity>16 & parity<666, 1, 0),
      agepreg1_outlier      =ifelse(age_preg1<12, 1,
                                    ifelse(age_preg1<666 & age_preg1>58, 1, 0)),
      agepreg2_outlier      =ifelse(age_preg2<12, 1,
                                    ifelse(age_preg2<666 & age_preg2>58, 1, 0)),
      agepreg3_outlier      =ifelse(age_preg3<12, 1,
                                    ifelse(age_preg3<666 & age_preg3>58, 1, 0)),
      agepreg4_outlier      =ifelse(age_preg4<12, 1,
                                    ifelse(age_preg4<666 & age_preg4>58, 1, 0)),
      agepreg5_outlier      =ifelse(age_preg5<12, 1,
                                    ifelse(age_preg5<666 & age_preg5>58, 1, 0)),
      agepreg6_outlier      =ifelse(age_preg6<12, 1,
                                    ifelse(age_preg6<666 & age_preg6>58, 1, 0)),
      agepreg7_outlier      =ifelse(age_preg7<12, 1,
                                    ifelse(age_preg7<666 & age_preg7>58, 1, 0)),
      agepreg8_outlier      =ifelse(age_preg8<12, 1,
                                    ifelse(age_preg8<666 & age_preg8>58, 1, 0)),
      agepreg9_outlier      =ifelse(age_preg9<12, 1,
                                    ifelse(age_preg9<666 & age_preg9>58, 1, 0)),
      agepreg10_outlier      =ifelse(age_preg10<12, 1,
                                     ifelse(age_preg10<666 & age_preg10>58, 1, 0)),
      breastfddur1_outlier  =ifelse(breastfeed_dur_b1<0,1,
                                    ifelse(breastfeed_dur_b1<666 & breastfeed_dur_b1>24, 1, 0)),
      breastfddur2_outlier  =ifelse(breastfeed_dur_b2<0,1,
                                    ifelse(breastfeed_dur_b2<666 & breastfeed_dur_b2>24, 1, 0)),
      breastfddur3_outlier  =ifelse(breastfeed_dur_b3<0,1,
                                    ifelse(breastfeed_dur_b3<666 & breastfeed_dur_b3>24, 1, 0)),
      breastfddur4_outlier  =ifelse(breastfeed_dur_b4<0,1,
                                    ifelse(breastfeed_dur_b4<666 & breastfeed_dur_b4>24, 1, 0)),
      breastfddur5_outlier  =ifelse(breastfeed_dur_b5<0,1,
                                    ifelse(breastfeed_dur_b5<666 & breastfeed_dur_b5>24, 1, 0)),
      breastfddur6_outlier  =ifelse(breastfeed_dur_b6<0,1,
                                    ifelse(breastfeed_dur_b6<666 & breastfeed_dur_b6>24, 1, 0)),
      breastfddur7_outlier  =ifelse(breastfeed_dur_b7<0,1,
                                    ifelse(breastfeed_dur_b7<666 & breastfeed_dur_b7>24, 1, 0)),
      breastfddur8_outlier  =ifelse(breastfeed_dur_b8<0,1,
                                    ifelse(breastfeed_dur_b8<666 & breastfeed_dur_b8>24, 1, 0)),
      breastfddur9_outlier  =ifelse(breastfeed_dur_b9<0,1,
                                    ifelse(breastfeed_dur_b9<666 & breastfeed_dur_b9>24, 1, 0)),
      breastfddur10_outlier  =ifelse(breastfeed_dur_b10<0,1,
                                     ifelse(breastfeed_dur_b10<666 & breastfeed_dur_b10>24, 1, 0)),
      
      #Set OC use outlier start age to min age at menarche
      #OC use duration outlier to 40 years of use
      
      ocusedur_outlier      =ifelse(ocuse_dur<0, 1,
                                    ifelse(ocuse_dur==666| ocuse_dur==777| ocuse_dur==888, 0,
                                           ifelse(ocuse_dur>480, 1, 0))), 
      ocusestart_outlier    =ifelse(ocuse_start<8, 1, 
                                    ifelse(ocuse_start<666 & ocuse_start>58, 1, 0)),
      ocusestop_outlier     =ifelse(ocuse_stop<8, 1, 
                                    ifelse(ocuse_stop<666 & ocuse_stop>58, 1, 0)),
      
      #Set duration of HRT use to 30 years - recommendation is 1-5 years...but many were outside this range
      hrtdur_outlier        =ifelse(hrt_dur<0, 1, 
                                    ifelse(hrt_dur>360 & hrt_dur<666, 1, 0)),
      hrtepdur_outlier      =ifelse(hrtep_dur<0, 1, 
                                    ifelse(hrtep_dur>360 & hrtep_dur<666,1, 0)),
      hrteonlydur_outlier   =ifelse(hrteonly_dur<0, 1,
                                    ifelse(hrteonly_dur>360 & hrteonly_dur<666, 1, 0)),
      
      #Change PA mets outlier to 100 - used as outlier in this analysis: https://pmc.ncbi.nlm.nih.gov/articles/PMC4451435/ 
      pamets_outlier        =ifelse(pa_mets<=0, 1, 
                                    ifelse(pa_mets>100 & pa_mets<666, 1, 0)),
      papct_outlier         =ifelse(pa_pct<0, 1, 
                                    ifelse(pa_pct<666 & pa_pct>100, 1, 0)),
      sizeprimary1_outlier  =ifelse(size_primary1<=0, 1,
                                    ifelse(size_primary1<666 & size_primary1>100, 1, 0)),
      sizeprimary2_outlier  =ifelse(size_primary2<=0, 1,
                                    ifelse(size_primary2<666 & size_primary2>100, 1, 0)),
      ki67primary1_outlier  =ifelse(ki67_primary1<0, 1,
                                    ifelse(ki67_primary1<666 & ki67_primary1>100, 1, 0)),
      ki67primary2_outlier  =ifelse(ki67_primary2<0, 1,
                                    ifelse(ki67_primary2<666 & ki67_primary2>100, 1, 0)))
  
 
  #####################################################
  #7. Add any needed variables for running validation.
  #   Change all "888", etc. to . 
  ####################################################
  
  #Change duration to 0 if never user
  data2<-data1 %>%          
    mutate(
      alcohol_dur=ifelse(alcohol_dur==777, 0, alcohol_dur),
      smoking_dur=ifelse(smoking_dur==777, 0, smoking_dur),
      breastfeed_dur=ifelse(breastfeed_dur==777, 0, breastfeed_dur),
      breastfeed_dur_b1=ifelse(breastfeed_dur_b1==777,0,breastfeed_dur),
      breastfeed_dur_b2=ifelse(breastfeed_dur_b2==777,0,breastfeed_dur),
      breastfeed_dur_b3=ifelse(breastfeed_dur_b3==777,0,breastfeed_dur),
      breastfeed_dur_b4=ifelse(breastfeed_dur_b4==777,0,breastfeed_dur),
      breastfeed_dur_b5=ifelse(breastfeed_dur_b5==777,0,breastfeed_dur),
      breastfeed_dur_b6=ifelse(breastfeed_dur_b6==777,0,breastfeed_dur),
      breastfeed_dur_b7=ifelse(breastfeed_dur_b7==777,0,breastfeed_dur),
      breastfeed_dur_b8=ifelse(breastfeed_dur_b8==777,0,breastfeed_dur),
      breastfeed_dur_b9=ifelse(breastfeed_dur_b9==777,0,breastfeed_dur),
      breastfeed_dur_b10=ifelse(breastfeed_dur_b10==777,0,breastfeed_dur),
      ocuse_dur=ifelse(ocuse_dur==777, 0, ocuse_dur),
      othcontracep_current=ifelse(othcontracep_current==777, 0, othcontracep_current),
      hrt_dur=ifelse(hrt_dur==777, 0, hrt_dur),
      hrtep_dur=ifelse(hrtep_dur==777,0, hrtep_dur),
      hrteonly_dur=ifelse(hrteonly_dur==777, 0, hrteonly_dur))
  
  #From here on can set all 666,777,888 to . for easier coding
  is.na(data2)<-data2==888
  is.na(data2)<-data2==777
  is.na(data2)<-data2==666
  
  #New dataset replacing outliers with NA for all variables - rename as data3 to preserve original
  data3<- data2 %>% 
    mutate( #outliers to NA
      waist=ifelse(waist_outlier==1, NA, waist),
      hip=ifelse(hip_outlier==1, NA, hip),
      age=ifelse(age_outlier ==1,NA, age),
      weight=ifelse(weight_outlier ==1, NA, weight),
      height=ifelse(height_outlier ==1 | height>200, NA, height),
      agemenarche=ifelse(agemenarche_outlier ==1, NA, agemenarche),
      bmi_earlyadult=ifelse(bmiearly_outlier ==1, NA, bmi_earlyadult),  
      bmi=ifelse(bmi_outlier ==1, NA, bmi),
      smoke_cigsperday=ifelse(smokingamt_outlier==1, NA, smoking_amt),
      agemenarche=ifelse(agemenarche_outlier==1, NA, agemenarche),
      parity=ifelse(parity_outlier==1, NA, parity),
      age_preg1=ifelse(agepreg1_outlier==1, NA, age_preg1),
      smoke_dur=ifelse(smokingdur_outlier==1, NA, smoking_dur),
      smoking_amt=ifelse(smokingamt_outlier==1, NA, smoking_amt),
      alc_dur=ifelse(alcoholdur_outlier==1, NA, alcohol_dur),
      alc_amt=ifelse(alcoholamt_outlier==1, NA, alcohol_amt),
      size_primary1=ifelse(sizeprimary1_outlier==1, NA, size_primary1),
      ki67_primary1=ifelse(ki67primary1_outlier==1, NA, ki67_primary1),
      ocuse_dur=ifelse(ocusedur_outlier==1, NA, ocuse_dur))
  
  
  #Create new variables for analysis and descriptive stats
  
  data4<-data3 %>%
    mutate(
      
      agemenarche_cat3=ifelse(agemenarche>0 & agemenarche<12, 1, 
                              ifelse(agemenarche >=12 & agemenarche<14, 2, 
                                     ifelse(agemenarche >=14, 3, NA))), 
      agemenarche_cat7=ifelse(agemenarche>0 & agemenarche<10.5, 1, 
                              ifelse(agemenarche>=10.5 & agemenarche<11.5, 2, 
                                     ifelse(agemenarche>=11.5 & agemenarche<12.5, 3, 
                                            ifelse(agemenarche>=12.5& agemenarche<13.5, 4, 
                                                   ifelse(agemenarche>=13.5&agemenarche<14.5, 5, 
                                                          ifelse(agemenarche>=14.5&agemenarche<15.5, 6, 
                                                                 ifelse(agemenarche>=15.5, 7, NA))))))),
      agemenarche_cat2=ifelse(agemenarche>0 & agemenarche<14, 1, 
                              ifelse(agemenarche>=14, 2, NA)),
      age_ge50=ifelse(age<50 & !is.na(age), 0, 
                      ifelse(age>=50, 1, NA)),
      age_cat1=ifelse(age>=45 & age<50, 1, 
                      ifelse(age>=50 & age<60,2,NA)),
      age_cat2=ifelse(age>=45 & age<50, 1, 
                      ifelse(age>=50&age<55, 2, 
                             ifelse(age>=55&age<60, 3, 
                                    ifelse(age>=60 & age<65, 4, 
                                           ifelse(age>=65 & age<75, 5, NA))))),
      height_in=height/2.54,
      height_cat3=ifelse(height>0 & height <160, 1,
                         ifelse(height >=160 & height<170, 2, 
                                ifelse(height>=170, 3, NA))),
      bmi_cat3=ifelse(bmi>0 & bmi<25, 1,
                      ifelse(bmi>=25 & bmi<30, 2, 
                             ifelse(bmi>=30, 3, NA))),
      bmi_cat4=ifelse(bmi>0 & bmi<18.5, 1, 
                      ifelse(bmi>=18.5&bmi<25, 2, 
                             ifelse(bmi>=25&bmi<30, 3, 
                                    ifelse(bmi>=30, 4, NA)))),
      bmi_cat5=ifelse(bmi>0 & bmi<21, 1, 
                      ifelse(bmi>=21&bmi<23, 2, 
                             ifelse(bmi>=23&bmi<35, 3, 
                                    ifelse(bmi>=25 & bmi<27, 4, 
                                           ifelse(bmi>=27, 5, NA))))),
      bmi_cat2=ifelse(bmi>0 & bmi<30, 1, 
                      ifelse(bmi>=30,2, NA )),
      
      #Fixes needed for CSDLH
      #Fix meno status if age at menopause is available.
      meno_status=ifelse(!is.na(meno_age) & age>meno_age, 1, meno_status), #setting to postmenopausal if record age is > meno_age
      #Change perimenopausal coded women to postmenopausal if age >65
      meno_status=ifelse(meno_status==3 & age>65, 1, meno_status),
      #code peri as premeno: 1 is now pre and 2 is post
      meno_status2=ifelse(meno_status==2 | meno_status==3, 1, ifelse(meno_status==1, 2, NA)),
      #use meno_status2 moving forward
      
      agemeno_combo=ifelse(meno_status2==1 & age_cat1==1, 1, 
                           ifelse(meno_status2==1 & age_cat1==2, 2, 
                                  ifelse(meno_status2==2 & age_cat2==1, 3, 
                                         ifelse(meno_status2==2 & age_cat2==2, 4, 
                                                ifelse(meno_status2==2 &age_cat2==3, 5, 
                                                       ifelse(meno_status2==2 &age_cat2==4, 6,NA)))))),
      
      ageflb=ifelse(age_preg1<=age & !is.na(age) & !is.na(age_preg1) , age_preg1, NA),
      
      #recode parous if parous and ageflb don't match - change parous to 1 if there is a valid age at first birth that fits in range above
      #recode parous if parous and parity do not match
      parous=ifelse(parous==1| (!is.na(ageflb) |(!is.na(parity)& parity>0)), 1, 0),
      parity_new=ifelse(parous==0, 0,
                        ifelse(parous==1 & parity>0 & !is.na(parity), parity, 
                               ifelse(parous==1 & (parity==0 | is.na(parity)), 1, NA))),
      #use parity_new variable moving forward for parity
      parity_cat=ifelse(parity_new<3, parity_new, 3),
      
      ageflb_cat5=ifelse(parous==0, 1, 
                         ifelse(ageflb>0 & ageflb<20, 2, 
                                ifelse(ageflb>=20 & ageflb<25, 3, 
                                       ifelse(ageflb>=25 & ageflb<30, 4, 
                                              ifelse(ageflb>=30, 5, NA))))),
      ageflb_cat4=ifelse(parous==0, 1, 
                         ifelse(ageflb>0 &ageflb<25, 2, 
                                ifelse(ageflb>=25&ageflb<30, 3, 
                                       ifelse(ageflb>=30, 4, NA)))),
      ageflb_parous=ifelse(ageflb>0 & ageflb<20, 1, 
                           ifelse(ageflb>=20 & ageflb<25, 2, 
                                  ifelse(ageflb>=25&ageflb<30, 3, 
                                         ifelse(ageflb>=30, 4, NA)))),
      
      famhx_first=ifelse(fhx_fdr_brca==0, 0, 
                         ifelse(fhx_fdr_brca==1, 1, NA)),
      biopsy_num=ifelse(!is.na(biopsies_number), biopsies_number, 
                        ifelse(is.na(biopsies_number) & !is.na(biopsies_yesno) & biopsies_yesno==1, 1, 
                               ifelse(is.na(biopsies_number) & !is.na(biopsies_yesno) & biopsies_yesno==0, 0, NA))),
      biopsy_num=suppressWarnings(as.numeric(biopsy_num)),
      
      #atypical hyperplasia, or biop with no atypia
      atyp_hyp=ifelse(biopsy_num>0 & (bbd_type1==3 | bbd_type2==3 | bbd_type3==3 | bbd_type4==3), 1, 
                      ifelse(biopsy_num>0 & (bbd_type1==1|bbd_type1==2|bbd_type2==1|bbd_type2==2|bbd_type3==1|bbd_type3==2), 0, NA)),
      biop_noatyp=ifelse(biopsy_num>0 & atyp_hyp==0, 1, 
                         ifelse(biopsy_num>0 & atyp_hyp==1, 0, NA)),
      
      biopsy_cat3=ifelse(biopsies_yesno==0 | bbd_type1==1, 1, 
                         ifelse(biopsies_yesno==1 & bbd_type1 !=2, 2, 
                                ifelse(bbd_type1==2, 3, NA))),
      
      biopsy_age50=ifelse(biopsy_num==0 & age_ge50==0, 1, 
                          ifelse(biopsy_num==1&age_ge50==0, 2, 
                                 ifelse(biopsy_num>1 & age_ge50==0, 3, 
                                        ifelse(biopsy_num==0 & age_ge50==1, 4, 
                                               ifelse(biopsy_num==1 & age_ge50==1, 5, 
                                                      ifelse(biopsy_num>1&age_ge50==1, 6, NA)))))),
      
      #Fix postmeno duration
      postmeno_dur_hold=ifelse (meno_status2==2, age-meno_age, 
                                ifelse(meno_status2==1, 0, NA)),
      
      postmeno_dur=ifelse(postmeno_dur_hold<0, 0, postmeno_dur_hold),
      
      hrt_use=ifelse(hrtuse==0|meno_status2==1, 1, 
                     ifelse(hrtuse==1 & meno_status2==2, 2, 
                            ifelse(hrtuse==2 & meno_status2==2, 3, 
                                   ifelse(hrtuse==3 & meno_status2==2, 4, NA)))),
      hrt_use_icare=ifelse(hrt_use==1, 0, 
                           ifelse(hrt_use==2, 1,
                                  ifelse(hrt_use==3, 2, 
                                         NA))),
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
      
      #reclassify alcohol status if "unkn. current or former and has alcohol_amt list for current use, change to current
      alcohol_status=ifelse(alcohol_amt>0 & !is.na(alcohol_amt) & alcohol_status==3, 1, alcohol_status),
      
      #Only assign an amount if they listed current drinker or unknown current/never, otherwise 0 if not current drinker, NA if unk. 
      alcgm=ifelse(alcohol_status==2 | alcohol_status==4, 0,
                   ifelse((alcohol_status==1|alcohol_status==3) & !is.na(alcohol_amt), alcohol_amt, NA)),
      alcgm_cat4=ifelse(alcohol_status==2|alcohol_status==4, 1, 
                        ifelse(alcgm>0 & alcgm<11, 2,
                               ifelse(alcgm>=11 & alcgm<22, 3,
                                      ifelse(alcgm>=22, 4, NA)))),
      alcgm_cat7=ifelse(alcgm_cat4==1, 1, 
                        ifelse(alcgm>0& alcgm<5, 2,
                               ifelse(alcgm>=5&alcgm<15,3,
                                      ifelse(alcgm>=15 & alcgm<25, 4,
                                             ifelse(alcgm>=25&alcgm<35, 5,
                                                    ifelse(alcgm>=35 & alcgm<45, 6,
                                                           ifelse(alcgm>=45, 7, NA))))))),
      #ocuse as never, past, current
      ocuse=ifelse(ocuse_ever==0,1,
                   ifelse(ocuse_ever==1 & is.na(ocuse_current), 2,
                          ifelse(ocuse_ever==1 & ocuse_current==0, 3,
                                 ifelse(ocuse_current==1, 4, NA)))),
      
      #meno categories age for iCARE
      meno_age_cat5=ifelse(is.na(meno_age), NA, ifelse(meno_age<40, 1, ifelse(meno_age>=40 & meno_age<45, 2, ifelse(meno_age>=45 & meno_age<50, 3,ifelse(meno_age>=50 & meno_age<55, 4, 5))))),
      #BC variables
      
      case_invasive=ifelse(!is.na(case) & case==0, 0,
                           ifelse(!is.na(case) & case==1 & invasive_primary1==1, 1, 
                                  ifelse(!is.na(case) & case==1 & invasive_primary1==2, 0, 
                                         ifelse(!is.na(case) & case==1 & is.na(invasive_primary1), NA, 0)))),
      
      #Additional case variables and time variables
      timediff_control=suppressWarnings(as.numeric(lastfup-record_year)),
      timediff_case=suppressWarnings(as.numeric(dxyear_primary1-record_year)),
      timediff=ifelse(case==1, timediff_case, 
                      ifelse(case==0, timediff_control, NA)),
      
      #Create indicators for case 5y, 6y, and within 10 years, inclusive
      #This starts from 
      case_5y=ifelse(case==1 & timediff_case<=5 & !is.na(timediff_case), 1, 0),
      case_5yinv=ifelse(case_5y==1 & case_invasive==1 & !is.na(case_5y) & !is.na(case_invasive), 1, 0),
      
      case_6y=ifelse(case==1 & timediff_case<=6 & !is.na(timediff_case), 1, 0),
      case_6yinv=ifelse(case_6y==1 & case_invasive==1 & !is.na(case_6y) & !is.na(case_invasive), 1, 0),
      
      case_10y=ifelse(case==1 & timediff_case<=10 & !is.na(timediff_case), 1, 0),
      case_10yinv=ifelse(case_10y==1 & case_invasive==1 & !is.na(case_5y) & !is.na(case_invasive),1,0),
      
      #round down dx age
      dxage=floor(as.integer(dxyear_primary1-birth_year)),
      dxage=ifelse(is.na(dxage), age+timediff_case, dxage),
      
      caseage_flag=ifelse(dxage<20| dxage>90, 1, 0),
      
      #follow up flag
      fup_flag=ifelse(timediff_control<0 | timediff_case<0, 1, 0))

####################################################
#8. Remove >75y and <20 for validation - cannot be calculated
####################################################

valdata<-data4[which (data4$age>=20 & data4$age<=75),]

#################################################################
#9.Add new ID variable that is simple numeric for each individual
#################################################################
rowvec<-c(1:nrow(valdata))
valdata$newid<-rowvec

#########################
# 10. SAVE DATASET
##########################

#PUT YOUR CORRECT FILEPATH

write.csv(valdata, "FILEPATH/bcrpp_valdata.csv")

#The dataset for validation is "valdata". Variables should be column names, Ids should be rows.
#I have found that R markdown will not always work if you "source" this file as opposed
#to saving the dataset and reading it in, which is why I suggest saving the dataset here for validation. 