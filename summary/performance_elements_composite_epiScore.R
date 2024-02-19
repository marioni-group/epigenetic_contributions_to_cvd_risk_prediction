library(survival)
library(survminer)
library(dplyr)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)
library(pROC)
library(ROCR)
library(ggplot2)

######## Yipeng's function to predict the outcome at onset using survival probabilities estimated by Breslow estimator #################

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(dataDF$tte))
  closest <- as.integer(uniqueTimes)
  thresholdIndex <- match(threshold, closest)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$tte, dataDF$event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
  onsetPredictions <- 1 - survivalPredictions
  
  # event should be 0 if tte is > 10
  dataDF$event <- sapply(1:nrow(dataDF), function(i) {
    if (dataDF$tte[[i]] > threshold) {
      0
    } else {
      dataDF$event[[i]]
    }
  })
  
  auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$event)
  prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$event)
  roc <- pROC::roc(response = dataDF$event, predictor = onsetPredictions)
  list(cumulativeBaseHaz = cumulativeBaseHaz, true_y = dataDF$event, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

# Rank Inverse Based Normalisation of the data
transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

x = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/assign_age_sex_30-70_assumptions/protein_EpiScore_2.csv')


# a selection of models that I want to study
riskFactorsOnlyCoxPH = coxph(Surv(tte, event) ~ scale(age) + factor(sex), x)
assignCoxPH = coxph(Surv(tte, event) ~ scale(age) + factor(sex) + scale(assign), x)
epiCoxPH = coxph(Surv(tte, event) ~ scale(age) + factor(sex) + scale(summed_scores), x)
proteinCoxPH = coxph(Surv(tte, event) ~ scale(age) + factor(sex) + scale(assign) + scale(summed_scores), x)
assignAlone = coxph(Surv(tte, event) ~ scale(assign), x)
assignEpiScore = coxph(Surv(tte, event) ~ scale(assign) + scale(summed_scores), x)
summary(assignEpiScore)

# models = list(r = riskFactorsOnlyCoxPH, a = assignCoxPH, e = epiCoxPH, f = proteinCoxPH, assignAlone = assignAlone, assignEpiScore = assignEpiScore)

# I chose only a couple of them
models = list(r = riskFactorsOnlyCoxPH, a = assignCoxPH, e = epiCoxPH, f = proteinCoxPH)

# quickly convert list to a dataframe
listToDataframe = function(list) {
  models = colnames(list)
  data = data.frame("Sensitivity" =  numeric(), "Specificity" = numeric(), "Model" = character())
  
  for (model in models) {
    tmp_data = data.frame("Sensitivity" =  list[['sensitivities', model]], "Specificity" = list[['specificities', model]], "Model" = model)
    data = rbind(data, tmp_data)
  }
  
  return(data)
}
threshold = 10

testResults <- lapply(models, function(m) {predictCoxPHOnset(x, m, threshold)})

rocs <- sapply(testResults, function(r) {r$roc})
rocs_df = listToDataframe(rocs)
rocs_df[which(rocs_df$Model =='r'), "Model"] = "age, sex"
rocs_df[which(rocs_df$Model =='e'), "Model"] = "age, sex, episcores"
rocs_df[which(rocs_df$Model =='a'), "Model"] = "age, sex, assign"
rocs_df[which(rocs_df$Model =='f'), "Model"] = "age, sex, assign, episcores"

cbPalette <- c("#000000", "blue", "red", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# draw a couple of ROC curves
rocs_df %>%
  ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
  geom_line() +
  theme_light() +
  scale_colour_manual(values=cbPalette) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2)) +
  xlab("False Positive Rate (1 - Specificity)") +
  ylab("True Postive Rate (Sensitivity)")

onsets = lapply(models, function(m) {predictCoxPHOnset(x, m, 10)})
roc.test(onsets$a$roc, onsets$f$roc)


metrics_proteins = data.frame(probability=0, tp_null=0, tp_full=0, tp_diff=0, 
                              fp_null=0, fp_full=0, fp_diff=0,
                              tn_null=0, tn_full=0, tn_diff=0,
                              fn_null=0, fn_full=0, fn_diff=0,
                              sensitivity_null=0, sensitivity_full=0,
                              specificity_null=0, specificity_full=0,
                              PPV_null=0, PPV_full=0, 
                              NPV_null=0, NPV_full=0)

for (i in seq(10, 90, by=5)) {
  # print(table(onsets$assign_protein$true_y, ifelse(onsets$assign_protein$onsetPredictions>i*0.1, 1, 0)))
  # print(table(testResults$f$true_y, ifelse(onsets$assign_protein$onsetPredictions>i*0.1, 1, 0)))
  assign_tp = sum(testResults$a$true_y == 1 & testResults$a$onsetPredictions>i*0.01)
  tp = sum(testResults$f$true_y == 1 & testResults$f$onsetPredictions>i*0.01)
  tp_diff = tp - assign_tp
  
  assign_fp = sum(testResults$a$true_y == 0 & testResults$a$onsetPredictions>i*0.01)
  fp = sum(testResults$f$true_y == 0 & testResults$f$onsetPredictions>i*0.01)
  fp_diff = fp - assign_fp
  
  assign_tn = sum(testResults$a$true_y == 0 & testResults$a$onsetPredictions<=i*0.01)
  tn = sum(testResults$f$true_y == 0 & testResults$f$onsetPredictions<=i*0.01)
  tn_diff = tn - assign_tn
  
  assign_fn = sum(testResults$a$true_y == 1 & testResults$a$onsetPredictions<=i*0.01)
  fn = sum(testResults$f$true_y == 1 & testResults$f$onsetPredictions<=i*0.01)
  fn_diff = fn - assign_fn
  
  sensitivity_null = assign_tp / (assign_tp + assign_fn)
  sensitivity_full = tp / (tp + fn)
  
  specificity_null = assign_tn / (assign_tn + assign_fp)
  specificity_full = tn / (tn + fp)
  #The positive predictive value is the probability that following a positive test result, that individual will truly have that specific disease.
  PPV_null = assign_tp / (assign_tp + assign_fp)
  PPV_full = tp / (tp + fp)
  
  NPV_null = assign_tn / (assign_tn + assign_fn)
  NPV_full = tn / (tn + fn)
  
  metrics_proteins[nrow(metrics_proteins) + 1,] = c(probability = i/100, 
                                                    tp_null=assign_tp, tp_full=tp, tp_diff=tp_diff, 
                                                    fp_null=assign_fp, fp_full=fp, fp_diff=fp_diff, 
                                                    tn_null=assign_tn, tn_full=tn, tn_diff=tn_diff, 
                                                    fn_null=assign_fn, fn_full=fn, fn_diff=fn_diff,
                                                    sensitivity_null=sensitivity_null, sensitivity_full=sensitivity_full,
                                                    specificity_null=specificity_null, specificity_full=specificity_full,
                                                    PPV_null=PPV_null, PPV_full=PPV_full,
                                                    NPV_null=NPV_null, NPV_full=NPV_full)
}

metrics_proteins = metrics_proteins[-1,]
View(metrics_proteins)

