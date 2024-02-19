### outlier removal - lifted from https://github.com/davebraze/FDB1/blob/master/R/outliers.R ###
library(kinship2)
library(coxme)
library(readxl)

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

troponin_levels = read.csv('/Volumes/marioni-lab/Generation_Scotland_data/Troponin/GS20K_Troponin_all.PHE', sep='\t');
troponin_levels = subset(troponin_levels, !(is.na(Troponin_T) & is.na(Troponin_I)))
# troponin I
count(subset(troponin_levels, !is.na(Troponin_I))) # 19130
dim(subset(troponin_levels, !is.na(Troponin_I)))
dim(subset(troponin_levels, !is.na(Troponin_T)))
count(subset(troponin_levels, !is.na(Troponin_I) & Troponin_I <1.2)) # 4815

troponin_levels$cTnI_corrected = troponin_levels$Troponin_I
# there are elements less than 1.2
troponin_levels[which(!is.na(troponin_levels$cTnI_corrected) & troponin_levels$cTnI_corrected < 1.2), "cTnI_corrected"] = 0.6

### body ###
body <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/clinical/body.csv")

# Age, sex
agesex = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/clinical/agesex.csv');

# Glucose, cholesterol, sodium, potassium, urea, creatinine
biochemistry = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/clinical/biochemistry.csv');

# Blood pressure - systolic, diasolic, heart rate + basic stats
bphr = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/clinical/BPHR.csv');

# Diabetes - type and date of diagnosis, unique
diabetes = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/clinical/diabetes.csv');

disease = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/PCQ/disease.csv')
disease = disease[1:15]

# Rheumatoid Arthritis
ra = read_excel('/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/RA.xlsx');

medicinesv5 = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/PCQ/medicationv5.csv');
medicinesv5 = medicinesv5[,c("ID", "bp")]
medicinesv5[which(is.na(medicinesv5$bp)), "bp"] = 0

medicinesv2 = read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/medicinesv2.csv')
bp_medicines = union(medicinesv2, medicinesv5)
# SMID 
smid = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/clinical/SIMD.csv');

smoking_PCQ = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/PCQ/smoking.csv');

kin_model <- readRDS("/Volumes/marioni-lab/Ola/Source/cox/kinship_matrix_using_fixed_2022-01-28-pedigree.rds")

# family history of stroke ?
all_merged = agesex %>%
  left_join(body[c(1,4)], c("id" = "id")) %>%
  left_join(troponin_levels, c("id" = "id")) %>%
  left_join(biochemistry, c("id" = "ID")) %>%
  left_join(diabetes, c("id" = "id")) %>%
  left_join(bphr, c("id" = "id")) %>%
  left_join(smoking_PCQ, c("id" = "id")) %>%
  left_join(smid, c("id" = "id")) %>%
  left_join(disease, c("id" = "ID")) %>%
  left_join(bp_medicines, c("id" = "ID"))


all_merged[which(is.na(all_merged$bp)), "bp"] = 0
all_merged$rank = all_merged$rank / 100
all_merged$diabetic = ifelse(!is.na(all_merged$tname), 1, 0)
all_merged$history_CVD = ifelse(all_merged$heart_disease_M | all_merged$heart_disease_F | all_merged$heart_disease_B | all_merged$heart_disease_S |
                                  all_merged$stroke_M | all_merged$stroke_F | all_merged$stroke_B | all_merged$stroke_S, 1, 0)
all_merged$sex.x <- factor(all_merged$sex.x, levels=c('F','M'), labels=c(0, 1))

cox = all_merged[c(
  "id",
  "Troponin_T",
  "Troponin_I",
  "cTnI_corrected",
  "rank",
  "history_CVD",
  "diabetic",
  "pack_years",
  "avg_sys",
  "Total_cholesterol",
  "HDL_cholesterol",
  "bp",
  "bmi"
)]

# [1]  0 14

names(cox) = c("id", "Troponin_T", "Troponin_I", "cTnI_corrected", "simd", "history_cvd", "diabetic",
               "pack_years", "sys_bp", "total_cholesterol", "HDL_cholesterol", "bp_medicines", "bmi")
dim(subset(cox, duplicated(cox$id)))
# [1]  0 14

### load in target file and merge ###
tar = read.csv("/Volumes/marioni-lab/Ola/Lab/Test_sets/gs20ktargets.tsv", sep='\t')
dim(tar)
# [1] 18414     6
tar2 = subset(tar, duplicated(tar$Sample_Sentrix_ID))
dim(tar2)
# [1] 0 6

cox <- merge(cox, tar, by.x="id", by.y="Sample_Name")
dim(subset(cox, duplicated(cox$id)))
# [1]  0 14

par(mfrow=c(3,5))
for(i in c(2:13)){
  title <- names(cox)[i]
  hist(cox[,i], breaks=100, main=title)
}

cox = cox[c(
  "id",
  "ID",
  "age",
  "sex",
  "Troponin_T",
  "Troponin_I",
  "cTnI_corrected",
  "pack_years",
  "sys_bp",
  "total_cholesterol",
  "HDL_cholesterol",
  "bmi",
  "Set"
)]
dim(subset(cox, duplicated(cox$ID)))
# [1]  0 13

#duplicates created here
names(cox)[2] = "Sample_Sentrix_ID"
dim(subset(cox, duplicated(cox$Sample_Sentrix_ID)))
# [1] 56 13


### recode outliers to NA ###
episcores = cox
episcores$Troponin_T <- outlierTrim(transform(episcores$Troponin_T))
episcores$Troponin_I <- outlierTrim(transform(episcores$Troponin_I))
episcores$cTnI_corrected <- outlierTrim(transform(episcores$cTnI_corrected))
episcores$pack_years <- outlierTrim(log(episcores$pack_years+1))
episcores$total_cholesterol <- outlierTrim(episcores$total_cholesterol)
episcores$HDL_cholesterol <- outlierTrim(episcores$HDL_cholesterol)
episcores$sys_bp <- outlierTrim(episcores$sys_bp)
episcores$bmi <- ifelse(episcores$bmi<18 | episcores$bmi>50, NA, log(episcores$bmi))

par(mfrow=c(2,4))
for(i in c(5:12)){
  title <- names(episcores)[i]
  hist(episcores[,i], breaks=100, main=title)
}


for(i in c(5:12)){
  title=names(episcores)[i]
  tmp = episcores[c(1, 2, 3, 4, 13, i)]
  tmp = na.omit(tmp)
  tmp[,6] = scale(resid(lmekin(tmp[,6] ~ tmp$age + factor(tmp$sex) +  I(tmp$age^2) + (1|tmp$id), varlist = kin_model*2)))
  hist(tmp[,6], breaks=100, main=title)
  write.csv(tmp, paste0('/Volumes/marioni-lab/Ola/Lab/EpiScores/pheno/', names(episcores[i]), ".csv"), row.names=F)
}


### make sure residuals aren't wildly different from original phenotypes ###
# for(i in c(4:13,15:22)){
#   print(names(tmp8)[i])
#   print(cor(tmp8[,i], tmp9[,i], use="pairwise.complete.obs"))
# }

/Users/shirin/Projects/R/general/troponin_episcores/episcores_using_biglasso/prepare-pheno-2022.R 
