library(dplyr)
library(readxl)
library(janitor)
box::use(../assign/assign_vars[...])

path = ''

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

troponin_levels = read.csv(paste0(path, 'Generation_Scotland_data/Troponin/GS20K_Troponin_all.PHE'), sep='\t');
troponin_levels = subset(troponin_levels, !(is.na(Troponin_T) & is.na(Troponin_I)))
# troponin I
count(subset(troponin_levels, !is.na(Troponin_I))) # 19130
count(subset(troponin_levels, !is.na(Troponin_I) & Troponin_I <1.2)) # 4815

troponin_levels$cTnI_corrected = troponin_levels$Troponin_I
# there are elements less than 1.2
troponin_levels[which(!is.na(troponin_levels$cTnI_corrected) & troponin_levels$cTnI_corrected < 1.2), "cTnI_corrected"] = 0.6

# Glucose, cholesterol, sodium, potassium, urea, creatinine
biochemistry = read.csv(paste0(path, 'Generation_Scotland_data/clinical/biochemistry.csv'));

# Blood pressure - systolic, diasolic, heart rate + basic stats
bphr = read.csv(paste0(path, 'Generation_Scotland_data/clinical/BPHR.csv'));

# Diabetes - type and date of diagnosis, unique
diabetes = read.csv(paste0(path, 'Generation_Scotland_data/clinical/diabetes.csv'));

disease = read.csv(paste0(path, 'Generation_Scotland_data_Sep2021/PCQ/disease.csv'))
disease = disease[1:15]

# Rheumatoid Arthritis
ra = read_excel(paste0(path, 'Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/RA.xlsx'));
ra = subset(ra, !duplicated(ra$id))
# SMID 
smid = read.csv(paste0(path, 'Generation_Scotland_data/clinical/SIMD.csv'));

smoking_PCQ = read.csv(paste0(path, 'Generation_Scotland_data_Sep2021/PCQ/smoking.csv'));
smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==1 & is.na(smoking_PCQ$cigs_day), smoking_PCQ$packs_day * 20, smoking_PCQ$cigs_day)
smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==2, smoking_PCQ$packs_day * 20, smoking_PCQ$cigs_day)
smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==3, round_half_up((smoking_PCQ$packs_day * 20) / 2.5), smoking_PCQ$cigarettes_per_day)

medicinesv5 = read.csv(paste0(path, 'Generation_Scotland_data_Sep2021/PCQ/medicationv5.csv'));
medicinesv5 = medicinesv5[,c("ID", "bp")]
medicinesv5[which(is.na(medicinesv5$bp)), "bp"] = 0

medicinesv2 = read.csv(paste0(path, 'Ola/Lab/ASSIGN/medicinesv2.csv'))
bp_medicines = union(medicinesv2, medicinesv5)

target = readRDS(paste0(path, "Ola/Lab/Test_sets/GS20k_Targets.rds"))

# family history of stroke ?
all_merged = target  %>%
  left_join(troponin_levels, c("Sample_Name" = "id")) %>%
  left_join(biochemistry, c("Sample_Name" = "ID")) %>%
  left_join(diabetes, c("Sample_Name" = "id")) %>%
  left_join(bphr, c("Sample_Name" = "id")) %>%
  left_join(smoking_PCQ, c("Sample_Name" = "id")) %>%
  left_join(smid, c("Sample_Name" = "id")) %>%
  left_join(disease, c("Sample_Name" = "ID")) %>%
  left_join(bp_medicines, c("Sample_Name" = "ID")) %>%
  left_join(ra, c("Sample_Name" = "id")) #%>%

all_merged[which(is.na(all_merged$bp)), "bp"] = 0
all_merged$rank = all_merged$rank / 100
all_merged$diabetic = ifelse(!is.na(all_merged$tname), 1, 0)
all_merged$ra = ifelse(!is.na(all_merged$section), 1, 0)
all_merged$history_CVD = ifelse(all_merged$heart_disease_M | all_merged$heart_disease_F | all_merged$heart_disease_B | all_merged$heart_disease_S |
                                  all_merged$stroke_M | all_merged$stroke_F | all_merged$stroke_B | all_merged$stroke_S, 1, 0)
all_merged$sex.x <- factor(all_merged$sex.x, 
                           levels=c('F','M'), 
                           labels=c(0, 1))

cox = all_merged[c(
  "Sample_Name",
  "sex.x",
  "age.x",
  "Troponin_T",
  "Troponin_I",
  "cTnI_corrected",
  "rank",
  "history_CVD",
  "diabetic",
  "cigs_day",
  "cigarettes_per_day",
  "pack_years",
  "avg_sys",
  "ra",
  "Total_cholesterol",
  "HDL_cholesterol",
  "bp",
  "Set"
)]


names(cox) = c("id", "sex", "age", "Troponin_T", "Troponin_I", "cTnI_corrected", "simd", "history_cvd", "diabetic",
               "cigs_day_raw", "cigs_day", "pack_years", "sys_bp", "ra", "total_cholesterol", "HDL_cholesterol", "bp_medicines", "set")

cox$assign = mapply(get_assign_score, cox$sex, cox$age, cox$ra, cox$total_cholesterol, cox$HDL_cholesterol, 
                    cox$sys_bp, cox$diabetic, cox$history_cvd, cox$cigs_day, cox$simd/10, cox$bp_medicines)

#cox$assign = mapply(get_assign_score, rep(0,19188), rep(63,19188), rep(0,19188), rep(6.2, 19188), rep(2.5, 19188), rep(134, 19188), rep(0, 19188), rep(0, 19188), rep(2, 19188), rep(3.693, 19188))

write.csv(cox, paste0(path, 'Ola/Lab/ASSIGN/GS_assign_data_for_summary.csv'), row.names = F)

