#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Ola
library("survival")
library("optparse")
library("stringr")
library("imputeTS")
library("jsonlite")
library("dplyr")
library("ggplot2")
library("hrbrthemes")
library("viridis")
library("pROC")
library("janitor")


transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

# Get arguments
######################################################

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Summary/SCORE2/cox_settings.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)
settings
# Import + prep data
######################################################

seed <- 42
set.seed(seed) # Set seed to ensure fold variation minimised 

df <- read.csv(paste0(settings$out_dir, "cox_covars_episcores_W1.csv"), check.names = FALSE) # 1996
# Remove people that have missing or strange time-to-event values (negative)
df = df[!is.na(df$tte) & df$tte>0, ]
df = na.omit(df)
table(df$event)
# 0    1 
# 3240  419 

if (settings$transform == "rank") {
  df$assign = transform(df$assign)
} else if (settings$transform == "log") {
  df$assign = log(df$assign + 1)
}

include_score2 = 1
if(include_score2 == 1) {
  score2 = read.csv('/Volumes/marioni-lab/Ola/Lab/SCORE2/score2_69.csv')
  score2 = score2[c('id', 'score2')]
  df = merge(df, score2, by="id")
  df$t_score2 = transform(df$score2)
}

df$cTnI_corrected = transform(df$cTnI_corrected)

x = df[c("Sentrix_ID", "age", "sex", "assign", "tte", "event", "cTnI_corrected", "score2", "t_score2")]

troponin = read.csv(settings$scores.troponin)
colnames(troponin)[2] = "troponin"
composite_109 = read.csv(settings$scores.composite_109)
colnames(composite_109)[2] = "composite_109"

include_troponin_composite = 0
if (include_troponin_composite == 1) {
  full = merge(x, troponin[-3], by.x = "Sentrix_ID", by.y = "sample")
} else {
  full = x
}

full = merge(full, composite_109[-3], by.x = "Sentrix_ID", by.y = "sample")
copy = full

threshold = 10

assignCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign, full)
score2CoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + t_score2, full)
proteinCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign + composite_109, full)
proteinCoxPHscore2 = coxph(Surv(tte, event) ~ age + factor(sex) + t_score2 + composite_109, full)
troponinCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign + transform(cTnI_corrected), full)
troponinCoxPHscore2 = coxph(Surv(tte, event) ~ age + factor(sex) + t_score2 + transform(cTnI_corrected), full)
fullCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign + composite_109 + transform(cTnI_corrected), full)
fullCoxPHscore2 = coxph(Surv(tte, event) ~ age + factor(sex) + t_score2 + composite_109 + transform(cTnI_corrected), full)

summary(assignCoxPH)
summary(score2CoxPH)
summary(proteinCoxPH)
summary(proteinCoxPHscore2)
summary(troponinCoxPH)
summary(troponinCoxPHscore2)
summary(fullCoxPH)
summary(fullCoxPHscore2)

models_assign = list(null = assignCoxPH, epiScore = proteinCoxPH, troponin = troponinCoxPH, full = fullCoxPH)
models_score2 = list(null = score2CoxPH, epiScore = proteinCoxPHscore2, troponin = troponinCoxPHscore2, full = fullCoxPHscore2)


# GENERATE INCREMENTAL R2
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
  concordance = summary(mod)$concordance[[1]]
  table=data.frame(cbind(beta,hr,se,z,lci,uci,p,concordance))
  return(table)
}

predictCoxPHOnset <- function(dataDF, coxPHModel) {
  
  uniqueTimes <- sort(unique(dataDF$tte))
  closest <- as.integer(uniqueTimes)
  thresholdIndex <- match(threshold, closest)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$tte, dataDF$event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))# check what this is doing
  onsetPredictions <- 1 - survivalPredictions
  
  cstat <- concordance(coxPHModel)$concordance
  
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
  
  list(cumulativeBaseHaz = cumulativeBaseHaz, 
       onsetPredictions = onsetPredictions, 
       auc = auc, prauc = prauc, roc = roc, 
       concordance = cstat)
}

listToDataframe = function(list) {
  models = colnames(list)
  data = data.frame("Sensitivity" =  numeric(), "Specificity" = numeric(), "Model" = character())
  
  for (model in models) {
    tmp_data = data.frame("Sensitivity" =  list[['sensitivities', model]], "Specificity" = list[['specificities', model]], "Model" = model)
    data = rbind(data, tmp_data)
  }
  
  return(data)
}

gather_model_statistics = function(df, mod) {
  statistics = extract_coxme_table(mod)
  predictCoxPHOnset = predictCoxPHOnset(df, mod)
  is_troponin = grepl("cTnI", rownames(extract_coxme_table(mod))[4])
  is_composite = grepl("composite", rownames(extract_coxme_table(mod))[4])
  
  formatted_list = data.frame(
    "age" = round_half_up(statistics[1,"hr"],2),
    "P_age" = statistics[1,"p"],
    "sex"= round_half_up(statistics[2,"hr"],2),
    "P_sex" = statistics[2,"p"], 
    "score" = round_half_up(statistics[3,"hr"],2),
    "P_score" = statistics[3,"p"],
    "CVD_episcore" = ifelse(is_composite, round_half_up(statistics[4,"hr"],2), NA),
    "P_CVD_episcore" = ifelse(is_composite, statistics[4,"p"], NA),
    "CVD_troponinscore" = ifelse(is_troponin, round_half_up(statistics[4,"hr"],2), round_half_up(statistics[5,"hr"],2)),
    "P_CVD_troponinscore" = ifelse(is_troponin, statistics[4,"p"], statistics[5,"p"]),
    "AUC" = round_half_up(predictCoxPHOnset$auc,3),
    "PRAUC" = round_half_up(predictCoxPHOnset$prauc,3),
    "Concordance" = round_half_up(statistics[1,"concordance"],3)
  )
  
  return(formatted_list)
}

threshold = 10

# count cases and controls
dataDF = x
dataDF$event <- sapply(1:nrow(dataDF), function(i) {
  if (dataDF$tte[[i]] > threshold) {
    0
  } else {
    dataDF$event[[i]]
  }
})
table(dataDF$event)

testResults <- lapply(models_assign, function(m) {predictCoxPHOnset(full, m)})
testResultsScore2 <- lapply(models_score2, function(m) {predictCoxPHOnset(full, m)})

statisticsAssign <- lapply(models_assign, function(m) {gather_model_statistics(full, m)})
statisticsAssign = bind_rows(statisticsAssign)
statisticsScore2 <- lapply(models_score2, function(m) {gather_model_statistics(full, m)})
statisticsScore2 = bind_rows(statisticsScore2)
write.csv(statisticsAssign, paste0(settings$out_dir, "statistics_ASSIGN.csv"), quote = F)
write.csv(statisticsScore2, paste0(settings$out_dir, "statistics_SCORE2.csv"), quote = F)

roc.test(testResults$null$roc, testResults$epiScore$roc)
roc.test(testResults$null$roc, testResults$troponin$roc)
roc.test(testResults$null$roc, testResults$full$roc)

rocs <- sapply(testResults, function(r) {r$roc})
rocs_df = listToDataframe(rocs)
rocs_df[which(rocs_df$Model =='r'), "Model"] = "age, sex"
rocs_df[which(rocs_df$Model =='e'), "Model"] = "age, sex, episcores"
rocs_df[which(rocs_df$Model =='a'), "Model"] = "age, sex, assign"
rocs_df[which(rocs_df$Model =='f'), "Model"] = "age, sex, assign, episcores"
# rocs_df[which(rocs_df$Model =='assignAlone'), "Model"] = "assign"
# rocs_df[which(rocs_df$Model =='assignEpiScore'), "Model"] = "assign, episcores"
cbPalette <- c("#000000", "blue", "red", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf(file=paste0(settings$out_dir, "AUC_", threshold, ".pdf"))

rocs_df %>%
  ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
  geom_line() +
  theme_light() +
  scale_colour_manual(values=cbPalette) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2)) +
  xlab("False Positive Rate (1 - Specificity)") +
  ylab("True Postive Rate (Sensitivity)")

dev.off()

aucs <- sapply(testResultsScore2, function(r) {r$auc})
praucs <- sapply(testResultsScore2, function(r) {r$prauc})
cstats <- sapply(testResultsScore2, function(r) {r$concordance})
metricsTable <- data.frame(AUC = aucs, PRAUC = praucs, C = cstats)
write.table(metricsTable, paste0(settings$out_dir, "metricsTable_", threshold, "_AUC.txt"), quote = F)

