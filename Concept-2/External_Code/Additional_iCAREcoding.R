##################################
#Additional iCARE variable coding
#For Kaitlyn when creating ref dataset.

#After running Step1.DataCleaning.R you will have the dataset valdata,
#saved as "FILEPATH/bcrpp_valdata.csv"
#Please read in this dataset

read.csv("FILEPATH/bcrpp_valdata.csv")

data_icare <- valdata %>%  
  #infinity if outside followup
  mutate(height=height/10, 
         age=floor(age),
         agestart=age,
         ageend=pmin(age+5, age+timediff),
         age_at_menarche=relevel(as.factor(agemenarche_cat7), ref="4"),
         parity=as.factor(parity_cat),
         oc_current=ocuse_current,
         oc_ever=ocuse_ever,
         bbd=bbd_history,
         age_first_birth=ifelse(parity_new>0, ageflb_parous, 
                                ifelse(parity_new==0 | parous==0, 1, NA)),
         age_first_birth=as.factor(age_first_birth),
         famhist=famhx_first,
         alcohol_intake=as.factor(alcgm_cat7),
         hrt=as.factor(hrt_use_icare),
         hrt_type=pmh_type_icare)

###########################
# Divide for lt50 and ge50
###########################

data_icare_ge50<-data_icare[which (data_icare$age>=50),]
data_icare_lt50<-data_icare[which (data_icare$age<50),]

## Fit names to iCARE objects
## bmi curc differs for lt50 and ge50
data_icare_lt50<- data_icare_lt50 %>% 
  mutate(bmi_curc=ifelse(bmi_cat4==1, 0, 
                         ifelse(bmi_cat4==2, 1, 
                                ifelse(bmi_cat4==3, 2, 
                                       ifelse(bmi_cat4==4,3, NA)))),
         bmi_curc=relevel(as.factor(bmi_curc), ref="1"))

