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
  if (length(unique(data_gail$Race))>1){
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
  upper.limit.absrisk.cat = observed.prob.cat + 1.96 * sd.absrisk.cat
  lower.limit.absrisk.cat = observed.prob.cat - 1.96 * sd.absrisk.cat
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


AUCcalc_ageinc= function(dataset1){
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



AUCcalc_ageadjust=function(dataset1){
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

    linear.predictor.age1=na.omit(as.numeric(dataset1$lp[dataset1$agegrp==1]))
    linear.predictor.age2=na.omit(as.numeric(dataset1$lp[dataset1$agegrp==2]))
    linear.predictor.age3=na.omit(as.numeric(dataset1$lp[dataset1$agegrp==3]))
    linear.predictor.age4=na.omit(as.numeric(dataset1$lp[dataset1$agegrp==4]))
    linear.predictor.age5=na.omit(as.numeric(dataset1$lp[dataset1$agegrp==5]))
    
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
     #Compute 95% confidence interval for the AUC estimate
     upper.limit.auc.age1 = auc.age1 + 1.96 * sqrt(var_auc.age1)
     lower.limit.auc.age1 = auc.age1 - 1.96 * sqrt(var_auc.age1)
     auc.age1<-cbind(auc.age1,upper.limit.auc.age1,lower.limit.auc.age1, var_auc.age1)
    
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
     #Compute 95% confidence interval for the AUC estimate
     upper.limit.auc.age2 = auc.age2 + 1.96 * sqrt(var_auc.age2)
     lower.limit.auc.age2 = auc.age2 - 1.96 * sqrt(var_auc.age2)
     auc.age2<-cbind(auc.age2,upper.limit.auc.age2,lower.limit.auc.age2, var_auc.age2)
     
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
     #Compute 95% confidence interval for the AUC estimate
     upper.limit.auc.age3 = auc.age3 + 1.96 * sqrt(var_auc.age3)
     lower.limit.auc.age3 = auc.age3 - 1.96 * sqrt(var_auc.age3)
     auc.age3<-cbind(auc.age3,upper.limit.auc.age3,lower.limit.auc.age3, var_auc.age3)
     
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
     #Compute 95% confidence interval for the AUC estimate
     upper.limit.auc.age4 = auc.age4 + 1.96 * sqrt(var_auc.age4)
     lower.limit.auc.age4 = auc.age4 - 1.96 * sqrt(var_auc.age4)
     auc.age4<-cbind(auc.age4,upper.limit.auc.age4,lower.limit.auc.age4, var_auc.age4)
     
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
     upper.limit.auc.age5 = auc.age5 + 1.96 * sqrt(var_auc.age5)
     lower.limit.auc.age5 = auc.age5 - 1.96 * sqrt(var_auc.age5)
     auc.age5<-cbind(auc.age5,upper.limit.auc.age5,lower.limit.auc.age5, var_auc.age5)

     #print all AUCs and variances by age groups
      colnames(auc.age1)<-c("auc","upper limit", "lower limit","variance")
      colnames(auc.age2)<-c("auc","upper limit", "lower limit","variance")
      colnames(auc.age3)<-c("auc","upper limit", "lower limit","variance")
      colnames(auc.age4)<-c("auc","upper limit", "lower limit","variance")
      colnames(auc.age5)<-c("auc","upper limit", "lower limit","variance")
      auc.age.adjust<-rbind(auc.age1, auc.age2, auc.age3, auc.age4, auc.age5)
      charq1<-as.character(q1)
      charq2<-as.character(q2)
      charq3<-as.character(q3)
      charq4<-as.character(q4)
      rownames(auc.age.adjust)<-cbind(paste0("<=", charq1), paste0(">",charq1, "-",charq2), paste0(">",charq2, "-", charq3), paste0(">",charq3,"-", charq4), paste0(">",charq4))
      return(auc.age.adjust)
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

returnrelrisk=function(dataset1){
observed.outcome=as.numeric(dataset1$event)
mean_cohort=mean(observed.outcome)
predicted.risk=as.numeric(dataset1$absrisk_predicted)
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

#break by deciles, training set 
deciles<-quantile(dataset1$lp, c(.10, .20, .30, .40, .50, .60, .70, .80, .90)) 
d1<-deciles[1]
d2<-deciles[2]
d3<-deciles[3]
d4<-deciles[4]
d5<-deciles[5]
d6<-deciles[6]
d7<-deciles[7]
d8<-deciles[8]
d9<-deciles[9]

dataset1$decile<-NA
dataset1$decile[(dataset1$lp < d1)]<-1
dataset1$decile[(dataset1$lp>=d1 & dataset1$lp<d2)]<-2
dataset1$decile[(dataset1$lp>=d2 & dataset1$lp<d3)]<-3
dataset1$decile[(dataset1$lp>=d3 & dataset1$lp<d4)]<-4
dataset1$decile[(dataset1$lp>=d4 & dataset1$lp<d5)]<-5
dataset1$decile[(dataset1$lp>=d5 & dataset1$lp<d6)]<-6
dataset1$decile[(dataset1$lp>=d6 & dataset1$lp<d7)]<-7
dataset1$decile[(dataset1$lp>=d7 & dataset1$lp<d8)]<-8
dataset1$decile[(dataset1$lp>=d8 & dataset1$lp<d9)]<-9
dataset1$decile[(dataset1$lp>=d9)]<-10

#E/O by risk deciles based on linear predictor
samp.size.cat = as.numeric(table(dataset1$decile))
#Observed and predicted probabilities at the categories
observed.prob.cat = aggregate(observed.outcome ~ dataset1$decile,
                              FUN = mean)$observed.outcome
predicted.prob.cat = aggregate(predicted.risk ~ dataset1$decile,
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
chisq.absrisk = sum(((observed.prob.cat - predicted.prob.cat)^2)/variance.absrisk.cat)
#sigma2_bar = mean(variance.absrisk.cat)
mu_bar = mean(observed.prob.cat)
#####
#Add relative risk 
#######
observed.RR.cat= observed.prob.cat/mu_bar
predicted.RR.cat=relrisk_gail
#Compute the variance-covariance matrix for the absolute risks
variance.mat.absrisk = diag(variance.absrisk.cat)

#Compute the derivative matrix for log relative risk
#error in derivate_diagonal
#variance.mat.risk

derivative_offdiagonal = -1/(10 * mu_bar)
derivative_diagonal = diag(1/observed.prob.cat + derivative_offdiagonal)
derivative_diagonal[lower.tri(derivative_diagonal)] = derivative_offdiagonal
derivative_diagonal[upper.tri(derivative_diagonal)] = derivative_offdiagonal
derivative.mat = derivative_diagonal[-10,]
varcov.logRR.cat = derivative_diagonal %*% variance.mat.absrisk %*% derivative_diagonal
sd.logRR.cat = sqrt(diag(varcov.logRR.cat))
sigmainv_logRR = solve(derivative.mat %*% variance.mat.absrisk %*% t(derivative.mat))
diff.logRR = (log(observed.RR.cat) - log(predicted.RR.cat))[-10]
chisq.logRR = as.numeric(diff.logRR %*% sigmainv_logRR %*% t(t(diff.logRR)))
#p-value for relative risk test
PVAL_logRR = 1 - pchisq(chisq.logRR, 10-1) 
#Compute 95% confidence interval for relative risk in each category
upper.limit.RR.cat = exp(log(observed.RR.cat) + 1.96 * sd.logRR.cat)
lower.limit.RR.cat = exp(log(observed.RR.cat) - 1.96 * sd.logRR.cat)
#Get values- E/O by cat- RR
EO_RR_bycat<-cbind(observed.RR.cat, lower.limit.RR.cat, upper.limit.RR.cat, predicted.RR.cat)
return(EO_RR_bycat)
}

plotcalibration_abs= function(dataset1){
  maxx=max(dataset1$predicted.prob.cat)*100+0.2
  maxy=max(dataset1$upper.limit.absrisk.cat)*100+0.2
  plot_test<-ggplot(data=dataset1,aes(x=predicted.prob.cat*100,y=observed.prob.cat*100)) +
    geom_point(color="firebrick",size=4) +
    geom_abline(intercept=0, slope=1,linetype="dotted") + geom_errorbar(aes(ymin=lower.limit.absrisk.cat*100, ymax=upper.limit.absrisk.cat*100), size=0.3) +
    labs(x='Expected risk (%)', y='Observed risk (%)', title="Absolute risk calibration")+ 
    scale_x_continuous(limits=c(0,maxx))  + scale_y_continuous(limits=c(0,maxy)) +
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
    scale_x_continuous(limits=c(0,5))  + scale_y_continuous(limits=c(0,5)) +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x = element_text(vjust = 0, size = 14),
          axis.title.y = element_text(vjust = 2, size = 14), 
          axis.title=element_text(size=16), legend.position="none") 
  return(plot_test)
}      

plotcalibration_rel= function(dataset1){
  maxx=max(dataset1$predicted.RR.cat)+0.2
  maxy=max(dataset1$upper.limit.RR.cat)+0.2
  plot_test<-ggplot(data=dataset1, aes(x=predicted.RR.cat, y=observed.RR.cat)) +
    geom_point(color="firebrick",size=4) +
    geom_abline(intercept=0, slope=1,linetype="dotted")  + geom_errorbar(aes(ymin=lower.limit.RR.cat, ymax=upper.limit.RR.cat)) +
    labs(x='Expected Relative Risk', y='Observed Relative Risk', title="Relative risk calibration")+ 
    scale_x_continuous(limits=c(0,max(maxx)))  + scale_y_continuous(limits=c(0,(maxy)))  +
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
    scale_x_continuous(limits=c(0,max(5)))  + scale_y_continuous(limits=c(0,(5)))  +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x = element_text(vjust = 0, size = 14),
          axis.title.y = element_text(vjust = 2, size = 14), 
          axis.title=element_text(size=16), legend.position="none")
  return(plot_test)
}
