library(dplyr)
episcores_cTnI_corrected = read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/EpiScores/cluster/test/episcores_cTnI_corrected.csv')
episcores_cTnI = read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/EpiScores/cluster/test/episcores_Troponin_I.csv')
episcores_cTnT = read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/EpiScores/cluster/test/episcores_Troponin_T.csv')

target = readRDS('/Volumes/marioni-lab/Ola/Lab/Test_sets/GS20k_Targets.rds')

# Performance Check
######################################################

# GENERATE INCREMENTAL R2

join <- left_join(episcores_cTnI_corrected, target, by = c("Sample_Sentrix_ID"="Sample_Sentrix_ID"))
null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
full <- summary(lm(pheno ~ age + sex + episcore, data=join))$r.squared
full_summary = summary(lm(pheno ~ age + sex + episcore, data=join))
print(full_summary)

print('Incremental R2')
print(round(100*(full - null), 3))

print('null R2')
print(null)

print('full R2')
print(full)

results <- data.frame("trait" = settings$feature, "null_R2" = null, "full_R2" = full, "incremental_R2" = round(100*(full - null), 3))
write.csv(results, file = paste0(output_dir, "results_", settings$feature, ".csv"), row.names = FALSE)