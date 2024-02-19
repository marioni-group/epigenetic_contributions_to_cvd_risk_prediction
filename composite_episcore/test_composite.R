#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Ola
library("survival")
library("optparse")
library("stringr")
library("jsonlite")
library("dplyr")
library("ggplot2")
library("hrbrthemes")
library("viridis")
library("pROC")


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

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/final/cox_settings.json'
if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)
settings
# Import + prep data
######################################################

set.seed(42) # Set seed to ensure fold variation minimised 
seed <- 42

# folds <- read.delim(opt$folds, header = TRUE, row.names = 2)
df <- read.csv(paste0(settings$out_dir, "cox_covars_episcores_W1.csv"), check.names = FALSE) # 3711
# Remove people that have missing or strange time-to-event values (negative)
df = df[!is.na(df$tte) & df$tte>0, ]
df = na.omit(df)

if (settings$transform == "rank") {
  df$assign = transform(df$assign)
} else if (settings$transform == "log") {
  df$assign = log(df$assign + 1)
}

rownames(df) = df$Sentrix_ID

x = subset(df, select = -c(Sentrix_ID, id, set, event, assign, tte))
x = as.matrix(x)

filename <- paste0(settings$out_dir, "elnet_coefficients_", settings$run, ".csv")
weights <-read.csv(filename, check.names = FALSE)
cat("\nHave imported weights for CpGs.\n")


# Filter test data
x_test <- x[,weights$Variable]
x_test <- t(x_test)
pred <- x_test * weights$Coefficient
sum <- as.data.frame(colSums(pred))
export_sum = data.frame("sample" = rownames(sum) , "summed_scores" = sum, "trait" = "EpiScore")
colnames(export_sum) = c("sample", "summed_scores", "trait")

# Output
######################################################

write.csv(export_sum, paste0(settings$out_dir, "export_sum", settings$run, ".csv"), row.names = F)
export_sum = read.csv(paste0(settings$out_dir, "export_sum", settings$run, ".csv"))
# Performance Check
######################################################

# x = df[c("age", "sex", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "tte", "event")]
x = df[c("age", "sex", "assign", "tte", "event", "cTnI_corrected")]
x$summed_scores = export_sum$summed_scores

write.csv(x, paste0(settings$out_dir, "protein_EpiScore", settings$run, ".csv"), row.names = F)

threshold = 10

assignCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign, x)
proteinCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign + summed_scores, x)
troponinCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign + transform(cTnI_corrected), x)
fullCoxPH = coxph(Surv(tte, event) ~ age + factor(sex) + assign + summed_scores + transform(cTnI_corrected), x)

summary(assignCoxPH)
summary(proteinCoxPH)
summary(troponinCoxPH)
summary(fullCoxPH)

null <- summary(assignCoxPH)$rsq
full <- summary(proteinCoxPH)$rsq

print(round(100*(full - null), 3))

# models = list(r = riskFactorsOnlyCoxPH, a = assignCoxPH, e = epiCoxPH, f = proteinCoxPH, assignAlone = assignAlone, assignEpiScore = assignEpiScore)
# models = list(r = riskFactorsOnlyCoxPH, a = assignCoxPH, e = epiCoxPH, f = proteinCoxPH)
models = list(null = assignCoxPH, 
              epi = proteinCoxPH,
              troponin = troponinCoxPH,
              full = fullCoxPH)

# GENERATE INCREMENTAL R2

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


testResults <- lapply(models, function(m) {predictCoxPHOnset(x, m)})

roc.test(testResults$null$roc, testResults$epi$roc)
# DeLong's test for two correlated ROC curves
# 
# data:  testResults$null$roc and testResults$full$roc
# Z = -1.4514, p-value = 0.1467
# alternative hypothesis: true difference in AUC is not equal to 0
# 95 percent confidence interval:
# -0.01368436  0.00204002
# sample estimates:
# AUC of roc1 AUC of roc2 
# 0.7185868   0.7244090 

rocs <- sapply(testResults, function(r) {r$roc})
rocs_df = listToDataframe(rocs)
rocs_df[which(rocs_df$Model =='null'), "Model"] = "age, sex, ASSIGN"
rocs_df[which(rocs_df$Model =='epi'), "Model"] = "age, sex, assign, episcores"
rocs_df[which(rocs_df$Model =='troponin'), "Model"] = "age, sex, assign, troponin"
rocs_df[which(rocs_df$Model =='full'), "Model"] = "age, sex, ASSIGN, CVD EpiScore, cTnI"
# rocs_df[which(rocs_df$Model =='assignAlone'), "Model"] = "assign"
# rocs_df[which(rocs_df$Model =='assignEpiScore'), "Model"] = "assign, episcores"
cbPalette <- c("#000000", "red", "red", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf(file=paste0(settings$out_dir, "AUC_red", threshold, ".pdf"))

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

# plot(1-rocs[['specificities', 'r']], rocs[['sensitivities', 'r']], type = "l", xlab="False Positives", ylab = "True Postives", col = "black")
# lines(1-rocs[['specificities', 'e']], rocs[['sensitivities', 'e']], type = "l", col = "red")
# lines(1-rocs[['specificities', 'a']], rocs[['sensitivities', 'a']], type = "l", col = "green")
# lines(1-rocs[['specificities', 'f']], rocs[['sensitivities', 'f']], type = "l", col = "blue")
# 
# legend(0.6, 0.2, legend=c("age, sex", "age, sex, episcores", "age, sex, assign", "age, sex, assign, episcores"),
#        col=c("black", "red", "green", "blue"), lty=1:1, cex=0.8)

aucs <- sapply(testResults, function(r) {r$auc})
praucs <- sapply(testResults, function(r) {r$prauc})
cstats <- sapply(testResults, function(r) {r$concordance})
metricsTable <- data.frame(AUC = aucs, PRAUC = praucs, C = cstats)

write.table(metricsTable, paste0(settings$out_dir, "metricsTable_", threshold, "_AUC.txt"), quote = F)

