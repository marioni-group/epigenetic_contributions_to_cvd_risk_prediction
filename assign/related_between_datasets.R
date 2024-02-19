library(dplyr)

# ped = read.csv('/Cluster_Filespace/Marioni_Group/Ola/Data/fixed-ped-2022-01-28-pedigree.csv')
# target_w4 = readRDS('/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-samplesheet_v3.rds')
# target_w1_w3 = readRDS('/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS10k_Targets.rds')

ped = read.csv('/Volumes/marioni-lab/Ola/Source/cox/fixed-ped-2022-01-28-pedigree.csv')
target_w4 = read.csv('/Volumes/marioni-lab/Ola/Lab/Test_sets/w4-target.csv')
target_w1_w3 = readRDS('/Volumes/marioni-lab/Ola/Lab/Test_sets/GS10k_Targets.rds')

# w4 = subset(target_w4, !startsWith(target_w4$Sample_Name, "S"))
w4 = target_w4
w4$Sample_Name = as.integer(w4$Sample_Name)
w4$Set = "W4"
w4 = w4 %>% left_join(ped, c("Sample_Name" = "volid"))
w4 = w4[c("Sample_Name", "Set", "famid")]

w3 = subset(target_w1_w3, target_w1_w3$Set %in% "W3")
w3 = w3 %>% left_join(ped, c("Sample_Name" = "volid"))
w3 = w3[c("Sample_Name", "Set", "famid")]

w1 = subset(target_w1_w3, target_w1_w3$Set %in% "W1")
w1 = w1 %>% left_join(ped, c("Sample_Name" = "volid"))
w1 = w1[c("Sample_Name", "Set", "famid")]

w1_w3 = union(w1, w3)
w3_w4 = union(w3, w4)
w1_w4 = union(w1, w4)

not_unique_w1_w3_w4 = subset(w1_w3, (w1_w3$famid %in% w4$famid)) # 5046
dim(not_unique_w1_w3_w4) 
unique_w1_w3_w4 = subset(w1_w3, !(w1_w3$famid %in% w4$famid)) # 4491
dim(unique_w1_w3_w4) 
# > 5046 + 4491
# [1] 9537

not_unique_w4_w1_w3 = subset(w4, (w4$famid %in% w1_w3$famid)) # 7941
dim(not_unique_w4_w1_w3) 
unique_w4_w1_w3 = subset(w4, !(w4$famid %in% w1_w3$famid)) # 936
dim(unique_w4_w1_w3) 
#7941+936 = 8877

not_unique_w1_and_w3_w4 = subset(w1, (w1$famid %in% w3_w4$famid)) # 1828
dim(not_unique_w1_and_w3_w4) 
unique_w1_and_w3_w4 = subset(w1, !(w1$famid %in% w3_w4$famid)) # 3259
dim(unique_w1_and_w3_w4)
# 1828 + 3259 = 5087

not_unique_w3_w4_and_w1 = subset(w3_w4, (w3_w4$famid %in% w1$famid)) #  3436
dim(not_unique_w3_w4_and_w1)  
unique_w3_w4_and_w1 = subset(w3_w4, !(w3_w4$famid %in% w1$famid)) # 9891
dim(unique_w3_w4_and_w1)
# 3436+9891 = 13327

not_unique_w1_w4_w3 = subset(w1_w4, (w1_w4$famid %in% w3$famid)) # 6206
dim(not_unique_w1_w4_w3) 
unique_w1_w4_w3 = subset(w1_w4, !(w1_w4$famid %in% w3$famid)) # 7758
dim(unique_w1_w4_w3)  
# 6206+7758=13964

not_unique_w3_w1_w4 = subset(w3, (w3$famid %in% w1_w4$famid)) # 3235
dim(not_unique_w3_w1_w4) 
unique_w3_w1_w4 = subset(w3, !(w3$famid %in% w1_w4$famid)) # 1215
dim(unique_w3_w1_w4) 

#write.csv(unique_w3_w4_and_w1, '/Volumes/marioni-lab/Ola/Lab/Test_sets/Unique_W3_W4_and_W1.csv', row.names = F)
