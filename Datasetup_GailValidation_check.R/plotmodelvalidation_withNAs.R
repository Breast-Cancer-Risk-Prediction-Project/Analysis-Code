plotModelValidationwithNAs = function(study.data, category.results, validation.results, 
                               dataset = "Example Dataset",
                               model.name = "Example Model",
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
                               risk.score.plot.percent.smooth = 50){
  
  followup = validation.results$Adjusted_Followup
  observed.outcome = validation.results$Subject_Specific_Observed_Outcome
  linear.predictor = validation.results$Subject_Specific_Risk_Score
  model.disease.incidence.rates = 
    validation.results$Population_Incidence_Rate
  timeframe = validation.results$Risk_Prediction_Interval
  ages = validation.results$Study_Incidence_Rate[,1]
  study_incidence = validation.results$Study_Incidence_Rate[,2]
  
  study.entry.age = study.data$study.entry.age
  sampling.weights = study.data$sampling.weights
  
  if(is.null(sampling.weights))
    sampling.weights = rep(1,dim(study.data)[1])
  
  freq = 1/sampling.weights
  
  freq.cases = freq[observed.outcome == 1]
  freq.controls = freq[observed.outcome == 0]
  
  linear.predictor.cases = 
    linear.predictor[observed.outcome == 1]
  
  linear.predictor.controls = 
    linear.predictor[observed.outcome == 0]
  
  observed.prob.cat = 
    as.numeric(as.character(Category.Results$Observed_Absolute_Risk))
  predicted.prob.cat = 
    as.numeric(as.character(Category.Results$Predicted_Absolute_Risk))
  
  observed.RR.cat = 
    as.numeric(as.character(Category.Results$Observed_Relative_Risk))
  predicted.RR.cat = 
    as.numeric(as.character(Category.Results$Predicted_Relative_Risk))
  
  upper.limit.absrisk.cat = 
    as.numeric(as.character(Category.Results$CI_Absolute_Risk_Upper))
  lower.limit.absrisk.cat = 
    as.numeric(as.character(Category.Results$CI_Absolute_Risk_Lower))
  
  upper.limit.RR.cat = 
    as.numeric(as.character(Category.Results$CI_Relative_Risk_Upper))
  lower.limit.RR.cat = 
    as.numeric(as.character(Category.Results$CI_Relative_Risk_Lower))
  
  chisq.absrisk = 
    as.numeric(validation.results$Hosmer_Lemeshow_Results[["statistic"]])
  number.of.percentiles = 
    as.numeric(validation.results$Hosmer_Lemeshow_Results[["parameter"]])
  PVAL_absrisk = 
    as.numeric(validation.results$Hosmer_Lemeshow_Results[["p.value"]])
  
  chisq.logRR = as.numeric(validation.results$RR_test_result[["statistic"]])
  PVAL_logRR = as.numeric(validation.results$RR_test_result[["p.value"]])
  
  auc = validation.results$AUC
  upper.limit.auc = validation.results$CI_AUC[2]
  lower.limit.auc = validation.results$CI_AUC[1]
  
  Overall.Expected.to.Observed.Ratio = 
    validation.results$Overall_Expected_to_Observed_Ratio
  CI.Overall.Expected.to.Observed.Ratio = 
    validation.results$CI_Overall_Expected_to_Observed_Ratio
  
  
  #adjust x and y axes to draw absolute risk calibration plot
  if(length(x.lim.absrisk) < 2){
    x.lim.absrisk = c(min(lower.limit.absrisk.cat, 
                          predicted.prob.cat)*100,
                      max(upper.limit.absrisk.cat, 
                          predicted.prob.cat)*100)
  }
  
  if(length(y.lim.absrisk) < 2){
    y.lim.absrisk = c(min(lower.limit.absrisk.cat, 
                          predicted.prob.cat)*100,
                      max(upper.limit.absrisk.cat, 
                          predicted.prob.cat)*100)
  }
  
  #adjust x and y axes to draw relative risk calibration plot
  if(length(x.lim.RR) < 2){
    x.lim.RR = c(min(lower.limit.RR.cat, predicted.RR.cat),
                 max(upper.limit.RR.cat, predicted.RR.cat))
  }
  
  if(length(y.lim.RR) < 2){
    y.lim.RR = c(min(lower.limit.RR.cat, predicted.RR.cat),
                 max(upper.limit.RR.cat, predicted.RR.cat))
  }
  
  m <- rbind(c(1,2),c(3,4),c(5,5))
  layout(m)
  
  oldpar =  par(mar = rep(4,4))
  plotCI(predicted.prob.cat*100, observed.prob.cat *100, ui =
           upper.limit.absrisk.cat * 100, 
         li = lower.limit.absrisk.cat * 100, 
         xlab = x.lab.absrisk, ylab = y.lab.absrisk, xlim = 
           x.lim.absrisk, 
         ylim = y.lim.absrisk, col = "red",
         pch = 16, pty = "s",cex.lab = 1.2)
  abline(0,1,lty=2,col="black")
  #mtext(paste("Hosmer-Lemeshow p-value:",PVAL),side = 3, line = 1)
  mtext("Absolute Risk Calibration",side = 3, line =
          1,font = 4)
  
  plotCI(predicted.RR.cat, observed.RR.cat, ui =
           upper.limit.RR.cat, li = lower.limit.RR.cat, xlab =
           x.lab.RR, ylab = y.lab.RR, xlim = x.lim.RR, ylim = y.lim.RR,
         col = "red",
         pch = 16, pty = "s",cex.lab = 1.2)
  abline(0,1,lty=2,col="black")
  #mtext(paste("Hosmer-Lemeshow p-value:",PVAL),side = 3, line = 1)
  mtext("Relative Risk Calibration",side = 3, line =
          1,font = 4)
  
  plot(density(linear.predictor.controls, 
               bw = risk.score.plot.bandwidth,
               kernel = risk.score.plot.kernel,
               weights = freq.controls/sum(freq.controls),
               n = (risk.score.plot.percent.smooth/100) * 
                 length(linear.predictor.controls)),
       xlim = c(min(density(linear.predictor.controls)$x,
                    density(linear.predictor.cases)$x),
                max(density(linear.predictor.controls)$x,
                    density(linear.predictor.cases)$x)),
       ylim = c(min(density(linear.predictor.controls)$y,
                    density(linear.predictor.cases)$y),
                max(density(linear.predictor.controls)$y,
                    density(linear.predictor.cases)$y)),
       main = "",xlab = "Risk Score", pty = "s",cex.lab = 1.2)
  lines(density(linear.predictor.cases,
                bw = risk.score.plot.bandwidth,
                kernel = risk.score.plot.kernel,
                weights = freq.cases/sum(freq.cases),
                n = (risk.score.plot.percent.smooth/100) * 
                  length(linear.predictor.controls)),col = "red")
  legend("topright",legend = c("Controls","Cases"),
         col = c("black","red"),
         lty=1,cex = 0.75,y.intersp = 2, xjust = 0.5, 
         yjust = 0.5, x.intersp = 3,adj = c(0.3,0.6))
  mtext("Discrimination",side = 3, line =
          1,font = 4)
  
  plot(model.disease.incidence.rates[,1],
       predict(loess(model.disease.incidence.rates[,2] 
                     ~ model.disease.incidence.rates[,1])),
       log = "y", type = "l",
       xlim = c(18,max(model.disease.incidence.rates[,1])),
       ylim = c(min(model.disease.incidence.rates[,2],study_incidence) + 
                  10^(-10),
                max(model.disease.incidence.rates[,2],study_incidence)),
       xlab = "Age (in years)",ylab = "Incidence Rates",pty = "s",
       cex.lab = 1.2)
  lines(ages,study_incidence,col = "red")
  legend("bottomright",legend = c("Population","Study"),
         col = c("black","red"),lty=1,
         cex = 0.9,y.intersp = 2, x.intersp = 3,adj = c(0.3,0.6))
  mtext("Incidence Rates",side = 3, line =
          1,font = 4)
  
  plot(c(0,2),c(0,2),xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',
       pty = "s")
  text(x = 1, y = 2, paste("Dataset: ",dataset))
  text(x = 1, y = 1.8, paste("Model Name: ",model.name))
  text(x = 1, y = 1.6, paste("Risk Prediction Interval:",timeframe))
  text(x = 1, y = 1.4, paste("Number of subjects (cases): "
                             ,length(observed.outcome),"(",
                             sum(observed.outcome),")"))
  text(x = 1, y = 1.2, paste("Follow-up time (years) [mean,range]: [",
                             round(mean(followup),3),", (",
                             round(range(followup)[1],3),",",
                             round(range(followup)[2],3),")"," ]"))
  text(x = 1, y = 1, paste("Baseline age (years) [mean,range]: [", 
                           round(mean(study.entry.age),3),", (",
                           round(range(study.entry.age)[1],3),",",
                           round(range(study.entry.age)[2],3),")"," ]"))
  text(x = 1, y = 0.8, paste("E/O [Estimate, 95% CI]:"
                             ,"[",round(Overall.Expected.to.Observed.Ratio,3),
                             ",","(",round(CI.Overall.Expected.to.Observed.Ratio[1],3),
                             ",",round(CI.Overall.Expected.to.Observed.Ratio[2],3),
                             ")","]"))
  text(x = 0.3, y = 0.5, "Absolute Risk Calibration", font = 2)
  text(x = 0.3, y = 0.3, paste("HL Test, df:",
                               round(chisq.absrisk,3),",",number.of.percentiles))
  #text(x = 0.1, y = 0.3, paste("HL Test df:",number.of.percentiles))
  text(x = 0.3, y = 0.1, paste("p-value:",
                               noquote(ifelse(PVAL_absrisk < 
                                                1e-100,"<1e-100",
                                              format.pval(PVAL_absrisk,digits=20, 
                                                          eps = 1e-100, 
                                                          scientific = TRUE)))))
  text(x = 1, y = 0.5, "Relative Risk Calibration", font = 2)
  text(x = 1, y = 0.3, paste("Test, df:",
                             round(chisq.logRR,3),",",number.of.percentiles - 1))
  #text(x = 0.5, y = 0.3, paste("Test df:",number.of.percentiles - 1))
  text(x = 1, y = 0.1, paste("p-value:",
                             noquote(ifelse(PVAL_logRR < 1e-100,"<1e-100",
                                            format.pval(PVAL_logRR,digits=20,
                                                        eps = 1e-100, 
                                                        scientific = 
                                                          TRUE)))))
  text(x = 1.7, y = 0.5, "Model Discrimination", font = 2)
  text(x = 1.7, y = 0.3, paste("AUC est:",round(auc,3)))
  text(x = 1.7, y = 0.1, 
       paste("95% CI: (",round(lower.limit.auc,3),",",
             round(upper.limit.auc,3),")"))
  
  par(oldpar)
  
  NULL
}
