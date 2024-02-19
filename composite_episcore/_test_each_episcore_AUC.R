library(survival)
library(survminer)
library(dplyr)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

extract_coxme_table <- function (mod){
  beta <- mod$coefficients #$fixed is not needed
  hr <-exp(beta)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  lci = exp(beta-1*1.96*se)
  uci = exp(beta+1*1.96*se)
  table=data.frame(cbind(beta,hr,se,z,lci,uci,p))
  return(table)
}

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(dataDF$tte))
  closest <- as.integer(uniqueTimes)
  thresholdIndex <- match(threshold, closest)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$tte, dataDF$event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))# check what this is doing
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
  list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

# Rank Inverse Based Normalisation of the data
transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

ped = "/Volumes/marioni-lab/Ola/Source/cox/kinship_matrix_using_fixed_2022-01-28-pedigree.rds"
kin_model = readRDS(ped)

cox = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/troponin/train/cox_covars.csv')
min_set = subset(cox, !is.na(tte) & tte>0)
min_set = subset(min_set, tte<20)
min_set = subset(min_set, !is.na(assign))

episcores_w1_w3 = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Protein_projections/EpiScore_projections_GS_9537.csv', check.names = FALSE)
episcores_w4 = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Protein_projections/EpiScore_projections_W4_8877_220221.csv', check.names = FALSE)
names(episcores_w4)[1] = 'ID'
identical(colnames(episcores_w1_w3), colnames(episcores_w4))
episcores = union(episcores_w1_w3, episcores_w4)
rownames(episcores) = episcores$ID
significant = read.csv('/Users/shirin/Desktop/hazard_ratios_protein_episcores_all_waves_only_significant.csv', check.names = FALSE)
episcores = episcores[,colnames(episcores) %in% significant$Protein]
episcores$ID = rownames(episcores)
x = merge(min_set, episcores, by.x="Sentrix_ID", by.y="ID")

# merged = subset(merged, set == "W1")

start = length(min_set) + 1
end = ncol(x)

settings = data.frame("output" = "/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/each_episcore_AUC/assign_transformed/")

AUC_results = x[1,start:end]
AUC_results = as.data.frame(t(AUC_results))
AUC_results$age_sex = 0
AUC_results$age_sex_assign = 0
AUC_results$age_sex_episcores = 0
AUC_results$age_sex_assign_episcores = 0
AUC_results = AUC_results[,c(-1)]

i = 1
for (epi in x[,start:end]) {
  
  protein = colnames(x)[i+start-1]
  
  riskFactorsOnlyCoxPH = coxph(Surv(tte, event) ~ age + sex, x)
  assignCoxPH = coxph(Surv(tte, event) ~ age + sex + transform(assign), x)
  epiCoxPH = coxph(Surv(x$tte, x$event) ~ x$age + x$sex + epi)
  proteinCoxPH = coxph(Surv(x$tte, x$event) ~ x$age + x$sex + epi + transform(x$assign))
  
  models = list(r = riskFactorsOnlyCoxPH, a = assignCoxPH, e = epiCoxPH, f = proteinCoxPH) 
  
  threshold = 10
  testResults <- lapply(models, function(m) {predictCoxPHOnset(x, m, threshold)})
  
  rocs <- sapply(testResults, function(r) {r$roc})
  
  pdf(file=paste0(settings$output, "AUC_", protein, "_threshold_", threshold, ".pdf"))
  
  plot(1-rocs[['specificities', 'r']], rocs[['sensitivities', 'r']], type = "l", xlab="False Positives", ylab = "True Postives", col = "black")
  lines(1-rocs[['specificities', 'e']], rocs[['sensitivities', 'e']], type = "l", col = "red")
  lines(1-rocs[['specificities', 'a']], rocs[['sensitivities', 'a']], type = "l", col = "green")
  lines(1-rocs[['specificities', 'f']], rocs[['sensitivities', 'f']], type = "l", col = "blue")
  
  legend(0.6, 0.2, legend=c("age, sex", paste0("age, sex, ", protein), "age, sex, assign", paste0("age, sex, assign, ", protein)),
         col=c("black", "red", "green", "blue"), lty=1:1, cex=0.8)
  dev.off()
  
  aucs <- sapply(testResults, function(r) {r$auc})
  praucs <- sapply(testResults, function(r) {r$prauc})
  metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)
  AUC_results[protein, "age_sex"] = aucs["r"]
  AUC_results[protein, "age_sex_assign"] = aucs["a"]
  AUC_results[protein, "age_sex_episcores"] = aucs["e"]
  AUC_results[protein, "age_sex_assign_episcores"] = aucs["f"]
  
  write.table(metricsTable, paste0(settings$output, "metricsTable_", protein, "_threshold_", threshold, "_AUC.txt"), quote = F)
  
  i = i+1
}


AUC_results$diff = AUC_results$age_sex_assign_episcores - AUC_results$age_sex_assign

write.csv(AUC_results, paste0(settings$output, "summary.csv"))

#######################################################

protein = "4337-49"
epi = x[,protein]

riskFactorsOnlyCoxPH = coxph(Surv(tte, event) ~ age + sex, x)
assignCoxPH = coxph(Surv(tte, event) ~ age + sex + transform(assign), x)
epiCoxPH = coxph(Surv(x$tte, x$event) ~ x$age + x$sex + epi)
proteinCoxPH = coxph(Surv(x$tte, x$event) ~ x$age + x$sex + epi + transform(x$assign))

models = list(r = riskFactorsOnlyCoxPH, a = assignCoxPH, e = epiCoxPH, f = proteinCoxPH) 


list = c(1:15)
settings = data.frame("output" = "/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/each_episcore_AUC/thresholds_for_4337-49/")

threshold_results = data.frame("age_sex" = 0, "age_sex_assign" = 0, "age_sex_episcores" = 0, "age_sex_assign_episcores" = 0, "threshold" = 0 )

for (threshold in list) {

  testResults <- lapply(models, function(m) {predictCoxPHOnset(x, m, threshold)})
  
  rocs <- sapply(testResults, function(r) {r$roc})
  
  pdf(file=paste0(settings$output, "AUC_", protein, "_threshold_", threshold, ".pdf"))
  
  plot(1-rocs[['specificities', 'r']], rocs[['sensitivities', 'r']], type = "l", xlab="False Positives", ylab = "True Postives", col = "black")
  lines(1-rocs[['specificities', 'e']], rocs[['sensitivities', 'e']], type = "l", col = "red")
  lines(1-rocs[['specificities', 'a']], rocs[['sensitivities', 'a']], type = "l", col = "green")
  lines(1-rocs[['specificities', 'f']], rocs[['sensitivities', 'f']], type = "l", col = "blue")
  
  legend(0.6, 0.2, legend=c("age, sex", paste0("age, sex, ", protein), "age, sex, assign", paste0("age, sex, assign, ", protein)),
         col=c("black", "red", "green", "blue"), lty=1:1, cex=0.8)
  dev.off()
  
  aucs <- sapply(testResults, function(r) {r$auc})
  praucs <- sapply(testResults, function(r) {r$prauc})
  metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)
  threshold_results = threshold_results %>% add_row("age_sex" = aucs["r"], "age_sex_assign" = aucs["a"], 
                                                    "age_sex_episcores" = aucs["e"], "age_sex_assign_episcores" = aucs["f"], "threshold" = threshold)

  write.table(metricsTable, paste0(settings$output, "metricsTable_", protein, "_threshold_", threshold, "_AUC.txt"), quote = F)
}

threshold_results = threshold_results[c(-1),]
threshold_results$diff = threshold_results$age_sex_assign_episcores - threshold_results$age_sex_assign
plot(threshold_results$threshold, threshold_results$diff, xlab = "Time threshold (y)", ylab="AUC(full model) - AUC(null model)", main = "4337-49" )

write.csv(threshold_results, paste0(settings$output, "summary.csv"), )




