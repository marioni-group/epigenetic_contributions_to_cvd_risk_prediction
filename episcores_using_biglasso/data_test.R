#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena and Ola

library("optparse")
library("stringr")
library("imputeTS")
library("jsonlite")
library("dplyr")
# Get arguments
######################################################

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Cluster_Filespace/Marioni_Group/Ola/Code/general/toolbox/troponin_episcores/episcores_using_biglasso/settings/cTnI_corrected_outliers.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)

# Import + prep data
######################################################

meth_data <- readRDS(settings$o_test_df)
sample_data <- read.csv(settings$o_pheno_test)
output_dir <- settings$o_dir

cat("\nHave imported testing methylation + sample data.\n")

# Make sure pheno and sample order matches
sample_data = subset(sample_data, !duplicated(sample_data$Sample_Sentrix_ID))
rownames(sample_data) = sample_data$Sample_Sentrix_ID
cat("Pheno after filtering:\n")
cat(paste("\t", dim(sample_data)))

meth_data = subset(meth_data, !duplicated(rownames(meth_data)))
meth_data <- meth_data[which(rownames(meth_data) %in% rownames(sample_data)),]
cat("Methylation after filtering:\n")
cat(paste("\t", dim(meth_data)))
meth_data <- meth_data[match(rownames(sample_data), rownames(meth_data)),]
cat("\nHave matched meth + sample rows.\n")

w <- paste0(settings$o_dir, "elnet_coefficients_", settings$feature, ".csv") 
weights <- read.csv(w, row.names = 1)
cat("\nHave imported weights for CpGs.\n")


# Only keep CpG data for selected features and prep x and y
intercept <- weights["Intercept","Coefficient"] # Extract intercept
features <- weights[2:nrow(weights),"Coefficient",drop=FALSE] # Keep only CpGs
features <- features[rownames(features) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data

# Keep just CpGs
x_test <- meth_data[,match(rownames(features), colnames(meth_data))]

# Match x and y
x_test <- t(x_test)
y_test <- sample_data[,settings$feature,drop=FALSE]
y_test <- y_test[colnames(x_test),,drop=FALSE]

cat("\nHave prepped + scaled per CpG as per user's choice.\n")

#Predictions
######################################################
pred <- x_test * features[,"Coefficient"]
pred_pp <- colSums(pred)
pred_pp <- pred_pp + intercept
pred_df <- data.frame("episcore" = pred_pp, "pheno" = y_test[[settings$feature]], "trait" = settings$feature)
pred_df <- data.frame("Sample_Sentrix_ID" = rownames(pred_df), pred_df)
cat("\nObtained predictions! Exporting...\n")

# Output
######################################################

# Output table
write.csv(pred_df, file = paste0(settings$o_dir, "episcores_", settings$feature, ".csv"), row.names = FALSE)

# Performance Check
######################################################

# GENERATE INCREMENTAL R2

join <- left_join(pred_df, sample_data, by = "Sample_Sentrix_ID")
null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
full <- summary(lm(pheno ~ age + sex + episcore, data=join))$r.squared

print('Incremental R2')
print(round(100*(full - null), 3))

print('null R2')
print(null)

print('full R2')
print(full)

results <- data.frame("trait" = settings$feature, "null_R2" = null, "full_R2" = full, "incremental_R2" = round(100*(full - null), 3))

summary(lm(pheno ~ age + sex + episcore, data=join))

write.csv(results, file = paste0(settings$o_dir, "results_", settings$feature, ".csv"), row.names = FALSE)
