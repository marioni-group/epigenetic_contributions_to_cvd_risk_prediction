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

# Rank Based Inverse Normalisation of the data
transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

serialize_coxph_summary <- function(mod, protein) {
  coefs = extract_coxme_table(mod)
  hr <-coefs[2,2]
  p <-coefs[2,7]
  lci <-coefs[2,5]
  uci <-coefs[2,6]
  assumptions = cox.zph(mod)
  local=assumptions$table["scale(transform(epi))", 'p']
  global=assumptions$table["GLOBAL", 'p']
  ret = data.frame(hr, p, lci, uci, local, global, protein)
  return(ret)
}
path = "/Volumes/marioni-lab/Ola/"
ped = paste0(path, "Source/cox/kinship_matrix_using_fixed_2022-01-28-pedigree.rds")
cox_path = paste0(path, 'Lab/Cox/basic_dataset_no_assumptions/cox_covars.csv')
episcores_w1_w3_path = paste0(path, 'Lab/EpiScores/Protein_projections/EpiScore_projections_GS_9537.csv')
episcores_w4_path = paste0(path, 'Lab/EpiScores/Protein_projections/EpiScore_projections_W4_8877_220221.csv')
cell_comp_path = paste0(path, "Lab/Test_sets/blood_cell_composition.tsv")

kin_model = readRDS(ped)
cox = read.csv(cox_path)
cox = select(cox, c("id", "Sentrix_ID", "age", "sex", "assign", "Troponin_I", 
                    "cTnI_corrected", "Troponin_T", "event", "tte"))

wbc = read.csv(cell_comp_path, sep='\t')
wbc = select(wbc, c("ID", "CD8T", "CD4T", "NK", "Bcell", "Gran"))

episcores_w1_w3 = read.csv(episcores_w1_w3_path, check.names = FALSE)
episcores_w4 = read.csv(episcores_w4_path, check.names = FALSE)
names(episcores_w4)[1] = 'ID'
identical(colnames(episcores_w1_w3), colnames(episcores_w4))
episcores = union(episcores_w1_w3, episcores_w4)

merged = merge(cox, episcores, by.x="Sentrix_ID", by.y="ID")
merged = subset(merged, !is.na(merged$cTnI_corrected))
merged = merge(merged, wbc, by.x="Sentrix_ID", by.y="ID")

start = length(cox) + 1
end = ncol(merged) - 5 # wbc counts

hazard_ratios = data.frame(hr=numeric(), p=numeric(), 
                           lci=numeric(), uci=numeric(), 
                           local=numeric(), global=numeric(),
                           protein=numeric())

concordance = data.frame(c_null=numeric(),
                         c_null_t=numeric(),
                         c_mod=numeric(), 
                         c_troponin=numeric(), 
                         diff=numeric(), protein=numeric())

hazard_ratios_troponin = hazard_ratios

i = 1

for (epi in merged[,start:end]) {
  #epi = merged[,start]
  protein = colnames(merged)[start-1+i]
  
  null = coxph(Surv(merged$tte, merged$event) ~
                  scale(transform(merged$assign)) +
                  scale(merged$CD4T) +
                  scale(merged$CD8T) +
                  scale(merged$NK) +
                  scale(merged$Bcell) +
                  scale(merged$Gran)
                )
  
  null_t = coxph(Surv(merged$tte, merged$event) ~
                   scale(transform(merged$assign)) +
                   scale(transform(merged$cTnI_corrected)) +
                   scale(merged$CD4T) +
                   scale(merged$CD8T) +
                   scale(merged$NK) +
                   scale(merged$Bcell) +
                   scale(merged$Gran)
                 )
  
  mod = coxph(Surv(merged$tte, merged$event) ~
                scale(transform(merged$assign)) +
                scale(transform(epi)) +
                scale(merged$CD4T) +
                scale(merged$CD8T) +
                scale(merged$NK) +
                scale(merged$Bcell) +
                scale(merged$Gran)
              )
  
  mod_troponin = coxph(Surv(merged$tte, merged$event) ~
                         scale(transform(merged$assign)) +
                         scale(transform(epi)) +
                         scale(transform(merged$cTnI_corrected)) +
                         scale(merged$CD4T) +
                         scale(merged$CD8T) +
                         scale(merged$NK) +
                         scale(merged$Bcell) +
                         scale(merged$Gran)
                       )
  
  episcore = serialize_coxph_summary(mod, protein)
  hazard_ratios = rbind(hazard_ratios, episcore)
  
  troponin = serialize_coxph_summary(mod_troponin, protein)
  hazard_ratios_troponin = rbind(hazard_ratios_troponin, troponin)
  
  cidx = concordance(null, null_t, mod, mod_troponin)
  
  c_null = cidx$concordance["null"]
  c_null_t = cidx$concordance["null_t"]
  c_mod = cidx$concordance["mod"]
  c_troponin = cidx$concordance["mod_troponin"]
  diff = c_troponin - c_mod
  
  concordance = rbind(concordance, data.frame(protein, c_null, c_null_t, c_mod, c_troponin, diff))
  
  i = i+1
} 

path = '/Volumes/marioni-lab/Ola/Lab/Cox/basic_dataset_no_assumptions/'

f1 = 'HRs_assign_epi_wbc.csv'
f2 = 'HRs_assign_epi_cTnI_corrected_wbc.csv'
f3 = 'concordance_cTnI_corrected_wbc.csv'

write.csv(hazard_ratios, paste0(path, f1), row.names = F)
write.csv(hazard_ratios_troponin, paste0(path, f2), row.names = F)
write.csv(concordance, paste0(path, f3), row.names = F)

hazard_ratios = subset(hazard_ratios, p<0.05)
hazard_ratios_troponin = subset(hazard_ratios_troponin, p<0.05)

f1 = 'HRs_assign_epi_significant_wbc.csv'
f2 = 'HRs_assign_epi_cTnI_corrected_significant_wbc.csv'

write.csv(hazard_ratios, paste0(path, f1), row.names = F)
write.csv(hazard_ratios_troponin, paste0(path, f2), row.names = F)
