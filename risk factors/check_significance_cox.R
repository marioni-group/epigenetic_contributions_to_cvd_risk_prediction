library(survival)
library(survminer)
library(dplyr)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

extract_coxme_table <- function (mod){
  beta <- mod$coefficients 
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

# Rank Inverse Based Normalisation of the data
transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

ped = "/Volumes/marioni-lab/Ola/Source/cox/kinship_matrix_using_fixed_2022-01-28-pedigree.rds"
kin_model = readRDS(ped)
cox = read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/Full_Dataset.csv')
cox = subset(cox, !is.na(tte) & tte>0)
cox = subset(cox, !is.na(assign) & !is.na(ID)) # 15664
cox = subset(cox, cox$age>=30 & cox$age<70) # 12971 items
cox$assign = outlierTrim(cox$assign, 3) 
cox = subset(cox, !is.na(assign)) #12790

cTnT = subset(cox, !is.na(Troponin_T))

mod_assign = coxph(Surv(cox$tte, cox$event) ~ cox$age + factor(cox$sex) + transform(cox$assign))

mod_assign = coxph(Surv(cox$tte, cox$event) ~ transform(cox$assign))
summary(mod_assign)
cox.zph(mod_assign)

mod_assign = coxme(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        (1|cox$id), varlist = kin_model*2)
summary(mod_assign)
####### cTnI #######

mod_troponins = coxph(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$cTnI)))
summary(mod_troponins)
cox.zph(mod_troponins)

mod_troponins = coxme(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(log(cox$cTnI)) +
                        (1|cox$id), varlist = kin_model*2)
summary(mod_troponins)

####### cTnI CORRECTED #######

mod_troponins = coxph(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$cTnI_corrected)))
summary(mod_troponins)
cox.zph(mod_troponins)

mod_troponins = coxme(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$cTnI_corrected)) +
                        (1|cox$id), varlist = kin_model*2)
summary(mod_troponins)

extract_coxme_table(mod_troponins)


####### cTnT #######

mod_troponins = coxph(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$Troponin_T))
                      )
summary(mod_troponins)
cox.zph(mod_troponins)

mod_troponins = coxme(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$Troponin_T)) +
                        (1|cox$id), varlist = kin_model*2)
summary(mod_troponins)
extract_coxme_table(mod_troponins)
####### cTnI + cTnT #######

mod_troponins = coxph(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$cTnI_corrected)) +
                        as.numeric(transform(cox$Troponin_T))
)
summary(mod_troponins)
cox.zph(mod_troponins)

mod_troponins = coxme(Surv(cox$tte, cox$event) ~
                        transform(cox$assign) +
                        as.numeric(transform(cox$cTnI_corrected)) +
                        as.numeric(transform(cox$Troponin_T)) +
                        (1|cox$id), varlist = kin_model*2)
summary(mod_troponins)
extract_coxme_table(mod_troponins)
####### EpiScore for cTnI (cTnT not valid) #######

hist(as.numeric(cox$episcore.cTnI_corrected), breaks = 20)
mod_troponins_epi = coxme(Surv(cox$tte, cox$event) ~
                            transform(cox$assign) +
                            transform(cox$episcore.cTnI_corrected) +
                            (1|cox$id), varlist = kin_model*2)
summary(mod_troponins_epi)

mod_troponins_epi = coxph(Surv(cox$tte, cox$event) ~
                            transform(cox$assign) +
                            cox$episcore.cTnI_corrected
                          )
summary(mod_troponins_epi)
cox.zph(mod_troponins_epi)

null <- summary(lm(cTnI_corrected ~ age + sex, data=cox))$r.squared
full <- summary(lm(cTnI_corrected ~ age + sex + episcore.cTnI_corrected, data=cox))$r.squared
full_summary = summary(lm(cTnI_corrected ~ age + sex + episcore.cTnI_corrected, data=cox))
print(full_summary)

print('Incremental R2')
print(round(100*(full - null), 3))

print('null R2')
print(null)

print('full R2')
print(full)


####### Epi_corrected_cTnI + corrected_cTnI + cTnT #######

mod_troponins_epi = coxme(Surv(cox$tte, cox$event) ~
                            transform(cox$assign) +
                            (as.numeric(log(cox$cTnI_corrected))) +
                            (as.numeric(log(cox$Troponin_T))) +
                            (cox$episcore.cTnI_corrected) +
                            (1|cox$id), varlist = kin_model*2)
summary(mod_troponins_epi)

mod_troponins_epi = coxph(Surv(cox$tte, cox$event) ~
                            transform(cox$assign) +
                            (as.numeric(log(cox$cTnI_corrected))) +
                            (as.numeric(log(cox$Troponin_T))) +
                            (cox$episcore.cTnI_corrected))
summary(mod_troponins_epi)
cox.zph(mod_troponins_epi)

####### EpiScores for Lifestyle factors #######
####### Pack years #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.pack_years))

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.pack_years) +
              (1|cox$id), varlist = kin_model*2)

summary(mod)

####### Total_cholesterol #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.total_cholesterol))

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.total_cholesterol) +
              (1|cox$id), varlist = kin_model*2)

summary(mod)

####### HDL cholesterol #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.HDL_cholesterol))

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.HDL_cholesterol) +
              (1|cox$id), varlist = kin_model*2)

summary(mod)

####### Sys_BP #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.sys_bp))

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.sys_bp) +
              (1|cox$id), varlist = kin_model*2)

summary(mod)

####### Sys_BP + HDL cholesterol #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.sys_bp) +
              scale(cox$episcore.HDL_cholesterol))

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(cox$episcore.sys_bp) +
              outlierTrim(cox$episcore.HDL_cholesterol) +
              (1|cox$id), varlist = kin_model*2)

summary(mod)


####### Lifestyle factors #######
hist(cox$total_cholesterol)

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(cox$total_cholesterol)
            )

summary(mod)
cox.zph(mod) ## Model invalid!

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(cox$total_cholesterol) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)

####### total cholesterol #######

hist(cox$total_cholesterol)

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1))
)

summary(mod)
cox.zph(mod) ## model invalid!

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1)) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)

####### HDL cholesterol #######

hist(cox$HDL_cholesterol)

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(cox$HDL_cholesterol)
)

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(cox$HDL_cholesterol) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)

####### Sys_BP #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
             scale(outlierTrim(cox$sys_bp))
)

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(outlierTrim(cox$sys_bp)) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)

####### Sys_BP #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(outlierTrim(cox$sys_bp))
)

summary(mod)
cox.zph(mod)

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              scale(outlierTrim(cox$sys_bp)) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)

####### pack_years + sys_bp #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1)) +
              outlierTrim(cox$total_cholesterol)
)

summary(mod)
cox.zph(mod) ## model invalid!

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1)) +
              outlierTrim(cox$total_cholesterol) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)


####### pack_years + total_cholesterol + HDL_cholesterol #######

mod = coxph(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1)) +
              outlierTrim(cox$total_cholesterol) +
              outlierTrim(cox$HDL_cholesterol)
              
)

summary(mod)
cox.zph(mod) ## model invalid!

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1)) +
              outlierTrim(cox$total_cholesterol) +
              outlierTrim(cox$HDL_cholesterol) +
              (1|cox$id), varlist = kin_model*2)
summary(mod)

####### pack_years + total_cholesterol + HDL_cholesterol #######

mod = coxph(Surv(tte, event) ~
              transform(assign) +
              outlierTrim(log(pack_years+1)) +
              outlierTrim(total_cholesterol) +
              outlierTrim(HDL_cholesterol) +
              outlierTrim(sys_bp), data = cox
            
)

summary(mod)
cox.zph(mod) ## model invalid!

mod = coxme(Surv(cox$tte, cox$event) ~
              transform(cox$assign) +
              outlierTrim(log(cox$pack_years+1)) +
              outlierTrim(cox$total_cholesterol) +
              outlierTrim(cox$HDL_cholesterol) +
              outlierTrim(cox$sys_bp) +
              (1|cox$id), varlist = kin_model*2) 
summary(mod)

tmp_set = cox[c("tte", "event", "assign", "pack_years", "total_cholesterol", "HDL_cholesterol", "sys_bp")]
tmp_set$pack_years = outlierTrim(log(cox$pack_years+1))
tmp_set$assign = transform(cox$assign)
tmp_set$total_cholesterol = outlierTrim(cox$total_cholesterol)
tmp_set$HDL_cholesterol = outlierTrim(cox$HDL_cholesterol)
tmp_set$sys_bp = outlierTrim(cox$sys_bp)
tmp_set = na.omit(tmp_set)

smokers = subset(tmp_set, pack_years>0)
non_smokers = subset(tmp_set, pack_years<=0)

mod = coxph(Surv(tte, event) ~
              assign +
              pack_years +
              total_cholesterol +
              HDL_cholesterol +
              sys_bp, data = smokers
            
)                           
                                 
summary(mod)
cox.zph(mod)                        
ggcoxfunctional(mod)

mod = coxph(Surv(tte, event) ~
              assign +
              pack_years +
              total_cholesterol +
              HDL_cholesterol +
              sys_bp, data = non_smokers
            
)      
summary(mod)
cox.zph(mod)                        
ggcoxfunctional(mod)

ggcoxfunctional(mod)
martingale = residuals(mod, type="martingale")

ggcoxdiagnostics(mod, type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(mod, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(mod, type = "martingale", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxzph(cox.zph(mod))
