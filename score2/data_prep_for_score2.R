library(janitor)
library(ggplot2)
library(ggpubr)

box::use(./score2[...])
#calculating SCORE2

#1 read in the assign dataset
assign = read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/GS_assign_data.csv')

#2 subset to variables used by score2, check if all variables are correctly coded
# birth date - month, year
# i am not quite sure what to do with age - maybe will need to round it to years
# sex - male, female
# systolic blood pressure - numeric
# total cholesterol - numeric
# hdl-cholesterol - numeric
# ldl-cholesterol - numeric - not needed
# current-smoker - yes / no
# diabetic ??
score2 = assign[c("id", "sex", "age", "diabetic", "sys_bp", "total_cholesterol", "HDL_cholesterol", "assign")]

smoking_PCQ = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/PCQ/smoking.csv');
table(as.factor(smoking_PCQ$status))
smoking_PCQ = subset(smoking_PCQ, smoking_PCQ$status != "")
table(as.factor(smoking_PCQ$status))
smoking_PCQ$smoker = ifelse(smoking_PCQ$status=="smoker", 1, 0)
smoking_PCQ = smoking_PCQ[c("id", "smoker")]

score2 = merge(score2, smoking_PCQ, by="id")


score2_69 = subset(score2, age < 70) # 16934

score2_69$score2 = mapply(get_score2_score, 
                          score2_69$sex, 
                          score2_69$age, 
                          score2_69$total_cholesterol, 
                          score2_69$HDL_cholesterol, 
                          score2_69$sys_bp,
                          score2_69$smoker, 
                          "low")

score2_40_69 = subset(score2_69, age >= 40 & age < 70) # 11348

ggplot(score2_69, aes(x=assign, y=score2)) +
  geom_point() +
  geom_density_2d() +
  geom_smooth(method=lm, colour='red') +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 6, label.x = 0) + 
  xlab("ASSIGN") +
  ylab("SCORE2")

ggplot(score2_40_69, aes(x=assign, y=score2)) +
  geom_point() +
  geom_density_2d() +
  geom_smooth(method=lm, colour='red') +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 6, label.x = 0) + 
  xlab("ASSIGN") +
  ylab("SCORE2")


# zapisac to wiekopojmne dzielo na klastrze w dwoch wersjach - Age range for SCORE2 i 30 - 69 (ass)
write.csv(score2_69, '/Volumes/marioni-lab/Ola/Lab/SCORE2/score2_69.csv', row.names = F)
write.csv(score2_40_69, '/Volumes/marioni-lab/Ola/Lab/SCORE2/score2_40_69.csv', row.names = F)

#age range for ASSIGN
#30-74

#age range for SCORE2
#40–69 years


#3 write a module calculating the score without adjustment for risk area
## PAMIETAC ZEBY TESTOWWAC KORELACJE POMIEDZY ASSIGNEM I SCORE2 TYLKO U LUDZI PONIZEJ 70 LAT (69 MAX)
#4 adjust for risk area

#5 test against the online tool