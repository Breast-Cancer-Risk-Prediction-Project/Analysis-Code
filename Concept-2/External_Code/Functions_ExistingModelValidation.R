#Functions for Existing model calibration analysis

summarizerisks<-function(data, var){
  data%>%
    summarise(n= sum(!is.na({{var}})),
              na= sum(is.na({{var}})),
              mean = mean({{var}}, na.rm = TRUE),
              sd = sd({{var}}, na.rm = TRUE),
              min=min({{var}}, na.rm=TRUE),
              q25 = quantile({{var}}, 0.25, na.rm = TRUE),
              median=median({{var}}, na.rm=TRUE),
              q75 = quantile({{var}}, 0.75, na.rm = TRUE),
              max= max({{var}}, na.rm=TRUE))
}


summarizerisks_race<-function(data, var){
  if (length(unique(data$Race))>1){
  data %>% group_by(Race) %>%
    summarise(n= sum(!is.na({{var}})),
              na= sum(is.na({{var}})),
              mean = mean({{var}}, na.rm = TRUE),
              sd = sd({{var}}, na.rm = TRUE),
              min=min({{var}}, na.rm=TRUE),
              q25 = quantile({{var}}, 0.25, na.rm = TRUE),
              median=median({{var}}, na.rm=TRUE),
              q75 = quantile({{var}}, 0.75, na.rm = TRUE),
              max= max({{var}}, na.rm=TRUE))
  }
}

##########################################
# Calibration functions
##########################################

#####
#O/E#
#####

obsexpected_overall= function(dataset1){
  #Now absolute risk calibration:
  observed.outcome=as.numeric(dataset1$observed.outcome)
  mean_cohort=mean(observed.outcome)
  predicted.risk=as.numeric(dataset1$absrisk)
  #this is length of dataset- # subjects
  length(predicted.risk)
  mean(predicted.risk)
  p_E = mean(predicted.risk,na.rm=TRUE)
  p_O = mean(observed.outcome,na.rm=TRUE)
  #variance of observed risk
  var_obsrisk = (p_O * (1 - p_O))/length(predicted.risk)
  #variance of log expected by observed
  var_log_exp_by_obs = var_obsrisk/p_O^2
  exp_by_obs = p_E/p_O
  lCI_exp_by_obs = exp(log(exp_by_obs) - 1.96 * sqrt(var_log_exp_by_obs))
  rCI_exp_by_obs = exp(log(exp_by_obs) + 1.96 * sqrt(var_log_exp_by_obs))
  exp_by_obs_overall<-cbind(exp_by_obs, lCI_exp_by_obs, rCI_exp_by_obs)
  return(exp_by_obs_overall)
}

obsexpected_bycat_abs<-function(dataset1){
  observed.outcome=as.numeric(dataset1$observed.outcome)
  mean_cohort=mean(observed.outcome)
  predicted.risk=as.numeric(dataset1$absrisk)
  #this is length of dataset- # subjects
  length(predicted.risk)
  mean(predicted.risk)
  p_E = mean(predicted.risk,na.rm=TRUE)
  p_O = mean(observed.outcome,na.rm=TRUE)
  #variance of observed risk
  var_obsrisk = (p_O * (1 - p_O))/length(predicted.risk)
  #variance of log expected by observed
  var_log_exp_by_obs = var_obsrisk/p_O^2
  exp_by_obs = p_E/p_O
  lCI_exp_by_obs = exp(log(exp_by_obs) - 1.96 * sqrt(var_log_exp_by_obs))
  rCI_exp_by_obs = exp(log(exp_by_obs) + 1.96 * sqrt(var_log_exp_by_obs))
  exp_by_obs_overall<-cbind(exp_by_obs, lCI_exp_by_obs, rCI_exp_by_obs)
  
  #break by deciles
  deciles<-quantile(dataset1$absrisk, c(.10, .20, .30, .40, .50, .60, .70, .80, .90)) 
  d1<-deciles[1]
  d2<-deciles[2]
  d3<-deciles[3]
  d4<-deciles[4]
  d5<-deciles[5]
  d6<-deciles[6]
  d7<-deciles[7]
  d8<-deciles[8]
  d9<-deciles[9]
  
  dataset1$decile_abs<-NA
  dataset1$decile_abs[(dataset1$absrisk <= d1)]<-1
  dataset1$decile_abs[(dataset1$absrisk>d1 & dataset1$absrisk<=d2)]<-2
  dataset1$decile_abs[(dataset1$absrisk>d2 & dataset1$absrisk<=d3)]<-3
  dataset1$decile_abs[(dataset1$absrisk>d3 & dataset1$absrisk<=d4)]<-4
  dataset1$decile_abs[(dataset1$absrisk>d4 & dataset1$absrisk<=d5)]<-5
  dataset1$decile_abs[(dataset1$absrisk>d5 & dataset1$absrisk<=d6)]<-6
  dataset1$decile_abs[(dataset1$absrisk>d6 & dataset1$absrisk<=d7)]<-7
  dataset1$decile_abs[(dataset1$absrisk>d7 & dataset1$absrisk<=d8)]<-8
  dataset1$decile_abs[(dataset1$absrisk>d8 & dataset1$absrisk<=d9)]<-9
  dataset1$decile_abs[(dataset1$absrisk>d9)]<-10
  
  #E/O by risk deciles based on linear predictor
  samp.size.cat = as.numeric(table(dataset1$decile_abs))
  #Observed and predicted probabilities at the categories
  observed.prob.cat = aggregate(observed.outcome ~ dataset1$decile_abs,
                                FUN = mean)$observed.outcome
  predicted.prob.cat = aggregate(predicted.risk ~ dataset1$decile_abs,
                                 FUN = mean)$predicted.risk
  #compute the variance in each category (observed category)
  variance.absrisk.cat = (observed.prob.cat * 
                            (1 - observed.prob.cat))/samp.size.cat
  #compute the standard deviation in each cateogry
  sd.absrisk.cat = sqrt(variance.absrisk.cat)
  #variance of log exp by observed
  #-------------------------------------------------------------------------------------------------- added #
  var_log_exp_by_obs.cat = variance.absrisk.cat/observed.prob.cat^2
  
  #E/O by cat
  exp_by_obs.cat=predicted.prob.cat/observed.prob.cat
  #Compute 95% confidence interval for E/O in each category
  upper.limit.exp_by_obs.cat = exp(log(exp_by_obs.cat) + 1.96 * var_log_exp_by_obs.cat)
  lower.limit.exp_by_obs.cat = exp(log(exp_by_obs.cat) - 1.96 * var_log_exp_by_obs.cat)
  #-------------------------------------------------------------------------------------------------- added #
  #Get values- E/O by cat:
  EO_bycat_test<-cbind(exp_by_obs.cat, lower.limit.exp_by_obs.cat, upper.limit.exp_by_obs.cat)
  #compute the Hosmer-Lemeshow statistic
  chisq.absrisk = sum(((observed.prob.cat - predicted.prob.cat)^2)
                      /variance.absrisk.cat)
  #sigma2_bar = mean(variance.absrisk.cat)
  ########
  #p-value for Hosmer Lemeshow test
  PVAL_absrisk = 1 - pchisq(chisq.absrisk, 10)
  #Compute 95% confidence interval for observed absolute risk in each 
  #category
  #add p value to end of EO_bycat_test
  EO_bycat_pval<-cbind(EO_bycat_test,PVAL_absrisk)
  return (EO_bycat_pval)
}

AUCcalc= function(dataset1){
  observed.outcome=as.numeric(dataset1$observed.outcome)
  linear.predictor = as.numeric(dataset1$lp)
  linear.predictor.cases = linear.predictor[observed.outcome == 1]
  linear.predictor.controls = linear.predictor[observed.outcome == 0]
  
  #matrix of indicators of risk score in cases bigger than risk score in controls
  indicator = vapply(linear.predictor.cases, function(x) x > 
                       linear.predictor.controls,
                     logical(length(linear.predictor.controls)))
  
  #estimate the AUC in full cohort setting
  auc = mean(indicator)
  auc
  #compute E_S_0[I(S_1 > S_0)]
  mean_S0_indicator = apply(indicator,2,mean)
  #compute E_S_1[I(S_1 > S_0)]
  mean_S1_indicator = apply(indicator,1,mean)
  #compute variance of AUC
  var_auc = var(mean_S0_indicator)/length(mean_S0_indicator) + 
    var(mean_S1_indicator)/length(mean_S1_indicator)
  #Compute 95% confidence interval for the AUC estimate
  upper.limit.auc = auc + 1.96 * sqrt(var_auc)
  lower.limit.auc = auc - 1.96 * sqrt(var_auc)
  auc_test<-cbind(auc,lower.limit.auc,upper.limit.auc)
  return(auc_test)
}

#AUC calc with age incorporated


AUCcalc_absrisk= function(dataset1){
  observed.outcome=as.numeric(dataset1$observed.outcome)
  absrisk.predictor = as.numeric(dataset1$absrisk)
  absrisk.predictor.cases = absrisk.predictor[observed.outcome == 1]
  absrisk.predictor.controls = absrisk.predictor[observed.outcome == 0]
  
  #matrix of indicators of risk score in cases bigger than risk score in controls
  indicator = vapply(absrisk.predictor.cases, function(x) x > 
                       absrisk.predictor.controls,
                     logical(length(absrisk.predictor.controls)))
  
  #estimate the AUC in full cohort setting
  auc = mean(indicator)
  auc
  #compute E_S_0[I(S_1 > S_0)]
  mean_S0_indicator = apply(indicator,2,mean)
  #compute E_S_1[I(S_1 > S_0)]
  mean_S1_indicator = apply(indicator,1,mean)
  #compute variance of AUC
  var_auc = var(mean_S0_indicator)/length(mean_S0_indicator) + 
    var(mean_S1_indicator)/length(mean_S1_indicator)
  #Compute 95% confidence interval for the AUC estimate
  upper.limit.auc = auc + 1.96 * sqrt(var_auc)
  lower.limit.auc = auc - 1.96 * sqrt(var_auc)
  auc_test<-cbind(auc,lower.limit.auc,upper.limit.auc)
  return(auc_test)
}


#Create function based on absolute risk with age adjustment- IPW of age groups
AUCcalc_absrisk_ageadjust=function(dataset1){
    #First: data by 5 age groups - based on age distribution of study
    quantiles_age<-quantile(dataset1$study.entry.age, c(.20, .40, .60, .80), na.rm=TRUE)   
    q1<-quantiles_age[1]
    q2<-quantiles_age[2]
    q3<-quantiles_age[3]
    q4<-quantiles_age[4]
    
    dataset1$agegrp<-NA
    dataset1$agegrp[(dataset1$study.entry.age <= q1)]<-1
    dataset1$agegrp[(dataset1$study.entry.age>q1 & dataset1$study.entry.age<=q2)]<-2
    dataset1$agegrp[(dataset1$study.entry.age>q2 & dataset1$study.entry.age<=q3)]<-3
    dataset1$agegrp[(dataset1$study.entry.age>q3 & dataset1$study.entry.age<=q4)]<-4
    dataset1$agegrp[(dataset1$study.entry.age>q4)]<-5

    linear.predictor.age1=na.omit(as.numeric(dataset1$absrisk[dataset1$agegrp==1]))
    linear.predictor.age2=na.omit(as.numeric(dataset1$absrisk[dataset1$agegrp==2]))
    linear.predictor.age3=na.omit(as.numeric(dataset1$absrisk[dataset1$agegrp==3]))
    linear.predictor.age4=na.omit(as.numeric(dataset1$absrisk[dataset1$agegrp==4]))
    linear.predictor.age5=na.omit(as.numeric(dataset1$absrisk[dataset1$agegrp==5]))
    
     #Compute AUC and 95% confidence interval for the AUC estimate for each age grp
    
     observed.outcome.age1=as.numeric(dataset1$observed.outcome[dataset1$agegrp==1])
     predicted.risk.age1=as.numeric(dataset1$absrisk[dataset1$agegrp==1])
        
     linear.predictor.cases.age1 = linear.predictor.age1[observed.outcome.age1 == 1 ]
     linear.predictor.controls.age1 = linear.predictor.age1[observed.outcome.age1 == 0 ]
     #matrix of indicators of risk score in cases bigger than risk score in controls
     indicator.age1 = vapply(linear.predictor.cases.age1, function(x) x > 
     linear.predictor.controls.age1,
     logical(length(linear.predictor.controls.age1)))
           
     #estimate the AUC in full cohort setting
     auc.age1 = mean(indicator.age1)
     #compute E_S_0[I(S_1 > S_0)]
     mean_S0_indicator.age1 = apply(indicator.age1,2,mean)
              
     #compute E_S_1[I(S_1 > S_0)]
     mean_S1_indicator.age1 = apply(indicator.age1,1,mean)
          
     #compute variance of AUC
     var_auc.age1 = var(mean_S0_indicator.age1)/length(mean_S0_indicator.age1) + 
     var(mean_S1_indicator.age1)/length(mean_S1_indicator.age1)
    
     auc.age1<-cbind(auc.age1, var_auc.age1)
    
     #Age 2
     observed.outcome.age2=as.numeric(dataset1$observed.outcome[dataset1$agegrp==2])
     predicted.risk.age2=as.numeric(dataset1$absrisk[dataset1$agegrp==1])
     
     linear.predictor.cases.age2 = linear.predictor.age2[observed.outcome.age2 == 1 ]
     linear.predictor.controls.age2 = linear.predictor.age2[observed.outcome.age2 == 0 ]
     #matrix of indicators of risk score in cases bigger than risk score in controls
     indicator.age2 = vapply(linear.predictor.cases.age2, function(x) x > 
                               linear.predictor.controls.age2,
                             logical(length(linear.predictor.controls.age2)))
     
     #estimate the AUC in full cohort setting
     auc.age2 = mean(indicator.age2)
     #compute E_S_0[I(S_1 > S_0)]
     mean_S0_indicator.age2 = apply(indicator.age2,2,mean)
     
     #compute E_S_1[I(S_1 > S_0)]
     mean_S1_indicator.age2 = apply(indicator.age2,1,mean)
     
     #compute variance of AUC
     var_auc.age2 = var(mean_S0_indicator.age2)/length(mean_S0_indicator.age2) + 
       var(mean_S1_indicator.age2)/length(mean_S1_indicator.age2)
    
     auc.age2<-cbind(auc.age2, var_auc.age2)
     
     #Age 3
     observed.outcome.age3=as.numeric(dataset1$observed.outcome[dataset1$agegrp==3])
     predicted.risk.age3=as.numeric(dataset1$absrisk[dataset1$agegrp==1])
     
     linear.predictor.cases.age3 = linear.predictor.age3[observed.outcome.age3 == 1 ]
     linear.predictor.controls.age3 = linear.predictor.age3[observed.outcome.age3 == 0 ]
     #matrix of indicators of risk score in cases bigger than risk score in controls
     indicator.age3 = vapply(linear.predictor.cases.age3, function(x) x > 
                               linear.predictor.controls.age3,
                             logical(length(linear.predictor.controls.age3)))
     
     #estimate the AUC in full cohort setting
     auc.age3 = mean(indicator.age3)
     #compute E_S_0[I(S_1 > S_0)]
     mean_S0_indicator.age3 = apply(indicator.age3,2,mean)
     
     #compute E_S_1[I(S_1 > S_0)]
     mean_S1_indicator.age3 = apply(indicator.age3,1,mean)
     
     #compute variance of AUC
     var_auc.age3 = var(mean_S0_indicator.age3)/length(mean_S0_indicator.age3) + 
       var(mean_S1_indicator.age3)/length(mean_S1_indicator.age3)
     
     auc.age3<-cbind(auc.age3,var_auc.age3)
     
     #Age 4
     observed.outcome.age4=as.numeric(dataset1$observed.outcome[dataset1$agegrp==4])
     predicted.risk.age4=as.numeric(dataset1$absrisk[dataset1$agegrp==1])
     
     linear.predictor.cases.age4 = linear.predictor.age4[observed.outcome.age4 == 1 ]
     linear.predictor.controls.age4 = linear.predictor.age4[observed.outcome.age4 == 0 ]
     #matrix of indicators of risk score in cases bigger than risk score in controls
     indicator.age4 = vapply(linear.predictor.cases.age4, function(x) x > 
                               linear.predictor.controls.age4,
                             logical(length(linear.predictor.controls.age4)))
     
     #estimate the AUC in full cohort setting
     auc.age4 = mean(indicator.age4)
     #compute E_S_0[I(S_1 > S_0)]
     mean_S0_indicator.age4 = apply(indicator.age4,2,mean)
     
     #compute E_S_1[I(S_1 > S_0)]
     mean_S1_indicator.age4 = apply(indicator.age4,1,mean)
     
     #compute variance of AUC
     var_auc.age4 = var(mean_S0_indicator.age4)/length(mean_S0_indicator.age4) + 
       var(mean_S1_indicator.age4)/length(mean_S1_indicator.age4)
     
     auc.age4<-cbind(auc.age4, var_auc.age4)
     
     #Age 5
     observed.outcome.age5=as.numeric(dataset1$observed.outcome[dataset1$agegrp==5])
     predicted.risk.age5=as.numeric(dataset1$absrisk[dataset1$agegrp==1])
     
     linear.predictor.cases.age5 = linear.predictor.age5[observed.outcome.age5 == 1 ]
     linear.predictor.controls.age5 = linear.predictor.age5[observed.outcome.age5 == 0 ]
     #matrix of indicators of risk score in cases bigger than risk score in controls
     indicator.age5 = vapply(linear.predictor.cases.age5, function(x) x > 
                               linear.predictor.controls.age5,
                             logical(length(linear.predictor.controls.age5)))
     
     #estimate the AUC in full cohort setting
     auc.age5 = mean(indicator.age5)
     #compute E_S_0[I(S_1 > S_0)]
     mean_S0_indicator.age5 = apply(indicator.age5,2,mean)
     
     #compute E_S_1[I(S_1 > S_0)]
     mean_S1_indicator.age5 = apply(indicator.age5,1,mean)
     
     #compute variance of AUC
     var_auc.age5 = var(mean_S0_indicator.age5)/length(mean_S0_indicator.age5) + 
       var(mean_S1_indicator.age5)/length(mean_S1_indicator.age5)
     #Compute 95% confidence interval for the AUC estimate
     auc.age5<-cbind(auc.age5,var_auc.age5)

     #IPW to get age-adjusted AUC
     
     #print all AUCs and variances by age groups
     colnames(auc.age1)<-c("auc1","variance1")
     colnames(auc.age2)<-c("auc2","variance2")
     colnames(auc.age3)<-c("auc3","variance3")
     colnames(auc.age4)<-c("auc4","variance4")
     colnames(auc.age5)<-c("auc5","variance5")
     
     auc.age.adjust<-cbind(auc.age1, auc.age2, auc.age3, auc.age4, auc.age5)
     
     auc.age.adjust<-as.data.frame(auc.age.adjust)
     auc.age.adjust$weight1<-1/auc.age.adjust$variance1
     auc.age.adjust$weight2<-1/auc.age.adjust$variance2
     auc.age.adjust$weight3<-1/auc.age.adjust$variance3
     auc.age.adjust$weight4<-1/auc.age.adjust$variance4
     auc.age.adjust$weight5<-1/auc.age.adjust$variance5
     
     #combined AUC, combined SE, combined confidence interval
     comb_auc=(auc.age.adjust$auc1*auc.age.adjust$weight1+ 
                 auc.age.adjust$auc2*auc.age.adjust$weight2+
                 auc.age.adjust$auc3*auc.age.adjust$weight3+
                 auc.age.adjust$auc4*auc.age.adjust$weight4+
                 auc.age.adjust$auc5*auc.age.adjust$weight5)
     comb_auc=comb_auc/(auc.age.adjust$weight1+auc.age.adjust$weight2+auc.age.adjust$weight3+auc.age.adjust$weight4+auc.age.adjust$weight5)
     comb_se=sqrt(1/(auc.age.adjust$weight1+auc.age.adjust$weight2+auc.age.adjust$weight3+auc.age.adjust$weight4+auc.age.adjust$weight5))
     comb_lci=comb_auc-(1.96*comb_se)
     comb_rci=comb_auc+(1.96*comb_se)
     
     
     auc.age.adjustresult<-cbind(comb_auc, comb_se, comb_lci, comb_rci)
     auc.age.adjustresult<-as.data.frame(auc.age.adjustresult)
     return(auc.age.adjustresult)
}


#Abs risk by categories - abs

returnabsrisk_abs<-function(dataset1){
  observed.outcome=as.numeric(dataset1$observed.outcome)
  mean_cohort=mean(observed.outcome)
  predicted.risk=as.numeric(dataset1$absrisk)
  #this is length of dataset- # subjects
  length(predicted.risk)
  mean(predicted.risk)
  p_E = mean(predicted.risk,na.rm=TRUE)
  p_O = mean(observed.outcome,na.rm=TRUE)
  #variance of observed risk
  var_obsrisk = (p_O * (1 - p_O))/length(predicted.risk)
  #variance of log expected by observed
  var_log_exp_by_obs = var_obsrisk/p_O^2
  exp_by_obs = p_E/p_O
  exp_by_obs
  #0.9036319
  lCI_exp_by_obs = exp(log(exp_by_obs) - 1.96 * sqrt(var_log_exp_by_obs))
  rCI_exp_by_obs = exp(log(exp_by_obs) + 1.96 * sqrt(var_log_exp_by_obs))
  exp_by_obs_overall<-cbind(exp_by_obs, lCI_exp_by_obs, rCI_exp_by_obs)
  #break by deciles of absolute risk
  deciles_abs<-quantile(dataset1$absrisk, c(.10, .20, .30, .40, .50, .60, .70, .80, .90), na.rm=TRUE) 
  d1<-deciles_abs[1]
  d2<-deciles_abs[2]
  d3<-deciles_abs[3]
  d4<-deciles_abs[4]
  d5<-deciles_abs[5]
  d6<-deciles_abs[6]
  d7<-deciles_abs[7]
  d8<-deciles_abs[8]
  d9<-deciles_abs[9]
  
  dataset1$decile_abs<-NA
  dataset1$decile_abs[(dataset1$absrisk <= d1)]<-1
  dataset1$decile_abs[(dataset1$absrisk>d1 & dataset1$absrisk<=d2)]<-2
  dataset1$decile_abs[(dataset1$absrisk>d2 & dataset1$absrisk<=d3)]<-3
  dataset1$decile_abs[(dataset1$absrisk>d3 & dataset1$absrisk<=d4)]<-4
  dataset1$decile_abs[(dataset1$absrisk>d4 & dataset1$absrisk<=d5)]<-5
  dataset1$decile_abs[(dataset1$absrisk>d5 & dataset1$absrisk<=d6)]<-6
  dataset1$decile_abs[(dataset1$absrisk>d6 & dataset1$absrisk<=d7)]<-7
  dataset1$decile_abs[(dataset1$absrisk>d7 & dataset1$absrisk<=d8)]<-8
  dataset1$decile_abs[(dataset1$absrisk>d8 & dataset1$absrisk<=d9)]<-9
  dataset1$decile_abs[(dataset1$absrisk>d9)]<-10
  #E/O by risk deciles based on linear predictor
  samp.size.cat = as.numeric(table(dataset1$decile_abs))
  samp.size.cat
  #Observed and predicted probabilities at the categories
  observed.prob.cat = aggregate(observed.outcome ~ dataset1$decile_abs,
                                FUN = mean)$observed.outcome
  predicted.prob.cat = aggregate(predicted.risk ~ dataset1$decile_abs,
                                 FUN = mean)$predicted.risk
  #compute the variance in each category (observed category)
  variance.absrisk.cat = (observed.prob.cat * 
                            (1 - observed.prob.cat))/samp.size.cat
  #compute the standard deviation in each cateogry
  sd.absrisk.cat = sqrt(variance.absrisk.cat)
  #variance of log exp by observed
  #-------------------------------------------------------------------------------------------------- added #
  var_log_exp_by_obs.cat = variance.absrisk.cat/observed.prob.cat^2
  #E/O by cat
  exp_by_obs.cat=predicted.prob.cat/observed.prob.cat
  #Compute 95% confidence interval for E/O in each category
  upper.limit.exp_by_obs.cat = exp(log(exp_by_obs.cat) + 1.96 * var_log_exp_by_obs.cat)
  lower.limit.exp_by_obs.cat = exp(log(exp_by_obs.cat) - 1.96 * var_log_exp_by_obs.cat)
  #-------------------------------------------------------------------------------------------------- added #
  #Get values- E/O by cat:
  EO_bycat_test<-cbind(exp_by_obs.cat, upper.limit.exp_by_obs.cat, lower.limit.exp_by_obs.cat)
  EO_bycat_test
  #compute the Hosmer-Lemeshow statistic
  chisq.absrisk = sum(((observed.prob.cat - predicted.prob.cat)^2)
                      /variance.absrisk.cat)
  ########
  #p-value for Hosmer Lemeshow test
  PVAL_absrisk = 1 - pchisq(chisq.absrisk, 10) 
  #Compute 95% confidence interval for observed absolute risk in each 
  #category
  upper.limit.absrisk.cat = observed.prob.cat + 1.96 * sd.absrisk.cat
  lower.limit.absrisk.cat = observed.prob.cat - 1.96 * sd.absrisk.cat
  #add p value to end of EO_bycat_test
  absrisk_bycat_abs<-cbind(observed.prob.cat, lower.limit.absrisk.cat, upper.limit.absrisk.cat, predicted.prob.cat)
  absrisk_bycat_abs<-as.data.frame(absrisk_bycat_abs)
  return(absrisk_bycat_abs)
}



plotcalibration_abs= function(dataset1){
  maxx=max(dataset1$predicted.prob.cat)*100+0.2
  maxy=max(dataset1$upper.limit.absrisk.cat)*100+0.2
  maxval=max(maxx,maxy)
  plot_test<-ggplot(data=dataset1,aes(x=predicted.prob.cat*100,y=observed.prob.cat*100)) +
    geom_point(color="firebrick",size=4) +
    geom_abline(intercept=0, slope=1,linetype="dotted") + geom_errorbar(aes(ymin=lower.limit.absrisk.cat*100, ymax=upper.limit.absrisk.cat*100), size=0.3) +
    labs(x='Expected risk (%)', y='Observed risk (%)', title="Absolute risk calibration")+ 
    scale_x_continuous(limits=c(-0.05,maxval))  + scale_y_continuous(limits=c(-0.05,maxval)) +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x = element_text(vjust = 0, size = 14),
          axis.title.y = element_text(vjust = 2, size = 14), 
          axis.title=element_text(size=16), legend.position="none") 
  return(plot_test)
}      

plotcalibration_abs2= function(dataset1){
  plot_test<-ggplot(data=dataset1,aes(x=predicted.prob.cat*100,y=observed.prob.cat*100)) +
    geom_point(color="firebrick",size=4) +
    geom_abline(intercept=0, slope=1,linetype="dotted") + geom_errorbar(aes(ymin=lower.limit.absrisk.cat*100, ymax=upper.limit.absrisk.cat*100), size=0.3) +
    labs(x='Expected risk (%)', y='Observed risk (%)', title="Absolute risk calibration")+ 
    scale_x_continuous(limits=c(-0.05,5))  + scale_y_continuous(limits=c(-0.05,5)) +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x = element_text(vjust = 0, size = 14),
          axis.title.y = element_text(vjust = 2, size = 14), 
          axis.title=element_text(size=16), legend.position="none") 
  return(plot_test)
}      

################
# relative risk
################

returnrelrisk=function(dataset1){
  #create categories of the risk score
  observed.outcome = dataset1$observed.outcome
  study.entry.age = dataset1$study.entry.age
  study.exit.age = dataset1$study.exit.age
  time.of.onset = dataset1$time.of.onset
  linear.predictor=dataset1$lp
  sampling.weights = dataset1$sampling.weights
  #observed.followup = study.exit.age - study.entry.age
  #followup = observed.followup
  predicted.risk=as.numeric(dataset1$absrisk)
  predicted.risk.interval = rep(5,length(observed.outcome))
  freq=1/sampling.weights
  freq.cases=freq[observed.outcome==1]
  freq.controls=freq[observed.outcome==0]
  linear.predictor.cases = linear.predictor[observed.outcome == 1]
  linear.predictor.controls = linear.predictor[observed.outcome == 0]
  linear.predictor.cat=quantcut(linear.predictor, q=10)
  number.of.percentiles=10
  #compute sample size at each category
  samp.size.cat = as.numeric(table(linear.predictor.cat))
  
  #Observed and predicted probabilities at the categories
  observed.prob.cat = aggregate(observed.outcome ~ linear.predictor.cat,
                                FUN = mean)$observed.outcome
  predicted.prob.cat = aggregate(predicted.risk ~ linear.predictor.cat,
                                 FUN = mean)$predicted.risk
  
  #compute the variance in each category
  variance.absrisk.cat = (observed.prob.cat * 
                            (1 - observed.prob.cat))/samp.size.cat
  
  #compute the standard deviation in each cateogry
  sd.absrisk.cat = sqrt(variance.absrisk.cat)
  
  #compute the Hosmer-Lemeshow statistic
  chisq.absrisk = sum(((observed.prob.cat - predicted.prob.cat)^2)
                      /variance.absrisk.cat)
  
  #sigma2_bar = mean(variance.absrisk.cat)
  mu_bar = mean(observed.prob.cat)
  
  #Producing relative risk estimates
  
  observed.RR.cat = observed.prob.cat/mu_bar
  predicted.RR.cat = predicted.prob.cat/mean(predicted.prob.cat)
  
  #Compute the variance-covariance matrix for the absolute risks
  variance.mat.absrisk = diag(variance.absrisk.cat)
  
  #Compute the derivative matrix for log relative risk
  derivative_offdiagonal = -1/(number.of.percentiles * mu_bar)
  derivative_diagonal = diag(1/observed.prob.cat + derivative_offdiagonal)
  
  derivative_diagonal[lower.tri(derivative_diagonal)] = derivative_offdiagonal
  derivative_diagonal[upper.tri(derivative_diagonal)] = derivative_offdiagonal
  
  derivative.mat = derivative_diagonal[-number.of.percentiles,]
  
  varcov.logRR.cat = derivative_diagonal %*% variance.mat.absrisk %*% 
    derivative_diagonal
  
  sd.logRR.cat = sqrt(diag(varcov.logRR.cat))
  
  sigmainv_logRR = solve(derivative.mat %*% variance.mat.absrisk %*% 
                           t(derivative.mat))
  
  diff.logRR = (log(observed.RR.cat) - 
                  log(predicted.RR.cat))[-10]
  
  chisq.logRR = as.numeric(diff.logRR %*% sigmainv_logRR %*% t(t(diff.logRR)))
  
  #p-value for relative risk test
  PVAL_logRR = 1 - pchisq(chisq.logRR, 10 - 1) 
  
  #Compute 95% confidence interval for relative risk in each category
  upper.limit.RR.cat = exp(log(observed.RR.cat) + 1.96 * sd.logRR.cat)
  lower.limit.RR.cat = exp(log(observed.RR.cat) - 1.96 * sd.logRR.cat)
  
  relriskcal<-cbind(observed.RR.cat, lower.limit.RR.cat, upper.limit.RR.cat, predicted.RR.cat)
  relriskcal<-as.data.frame(relriskcal)
  return(relriskcal)
}



plotcalibration_rel= function(dataset1){
  maxx=max(dataset1$predicted.RR.cat)+0.2
  maxy=max(dataset1$upper.limit.RR.cat)+0.2
  plot_test<-ggplot(data=dataset1, aes(x=predicted.RR.cat, y=observed.RR.cat)) +
    geom_point(color="firebrick",size=4) +
    geom_abline(intercept=0, slope=1,linetype="dotted")  + geom_errorbar(aes(ymin=lower.limit.RR.cat, ymax=upper.limit.RR.cat)) +
    labs(x='Expected Relative Risk', y='Observed Relative Risk', title="Relative risk calibration")+ 
    scale_x_continuous(limits=c(-0.05,max(maxx)))  + scale_y_continuous(limits=c(-0.05,(maxy)))  +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x = element_text(vjust = 0, size = 14),
          axis.title.y = element_text(vjust = 2, size = 14), 
          axis.title=element_text(size=16), legend.position="none")
  return(plot_test)
}

plotcalibration_rel2= function(dataset1){
  plot_test<-ggplot(data=dataset1, aes(x=predicted.RR.cat, y=observed.RR.cat)) +
    geom_point(color="firebrick",size=4) +
    geom_abline(intercept=0, slope=1,linetype="dotted")  + geom_errorbar(aes(ymin=lower.limit.RR.cat, ymax=upper.limit.RR.cat)) +
    labs(x='Expected Relative Risk', y='Observed Relative Risk', title="Relative risk calibration")+ 
    scale_x_continuous(limits=c(-0.05,max(5)))  + scale_y_continuous(limits=c(-0.05,(5)))  +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x = element_text(vjust = 0, size = 14),
          axis.title.y = element_text(vjust = 2, size = 14), 
          axis.title=element_text(size=16), legend.position="none")
  return(plot_test)
}
