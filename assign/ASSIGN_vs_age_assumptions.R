library(dplyr)
library(optparse)
library(jsonlite)
library(survival)
library(survminer)
library(coxme)

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

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

make_cox_dataset = function(ds) {
  
  ds[which(is.na(ds$mob)), "mob"] = 1
  # there is another file - targets file - where this has been corrected
  
  ds$dead = 0
  # add natural deaths
  ds[which(!is.na(ds$dod_ym)), "dead"] = 1
  
  #check for multiple events!!!
  ds$event = 0
  ds[which(!is.na(ds$Incident.Event)), "event"] = 1
  
  # year and month of event
  ds$yoe = as.numeric(substring(ds$Event_Date.Event, 1, 4))
  ds$moe = as.numeric(substring(ds$Event_Date.Event, 5, 6))
  
  # year and month of death
  ds$yod = as.numeric(substring(ds$dod_ym, 1, 4))
  ds$mod = as.numeric(substring(ds$dod_ym, 5, 6))
  
  ds$censor_y = ifelse(ds$event == 1, ds$yoe, ifelse(ds$dead == 1, ds$yod, 2021))
  ds$censor_m = ifelse(ds$event == 1, ds$moe, ifelse(ds$dead == 1, ds$mod, 9))
  
  #censoring
  ds$yr_diff = ifelse(ds$event==0, ds$censor_y - ds$yob, ds$yoe - ds$yob)
  ds$m_diff = ifelse(ds$event==0, ((ds$censor_m - ds$mob)/12), ((ds$moe - ds$mob)/12))
  
  ds$age_event = ds$yr_diff + ds$m_diff
  ds$tte = ds$age_event - ds$age
  ds$tte = ifelse(ds$tte < -1, NA, ds$tte)
  ds$tte = ifelse(ds$tte < 0, 0, ds$tte)
  
  return(ds)
}

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/main_cox_model/cox_settings.json'
path = ''
if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)
settings

seed <- 42
set.seed(seed) # Set seed to ensure fold variation minimised 


####################################### Input data #######################################

agesex = read.csv(settings$prep.agesex);
agemonths = read.csv(settings$prep.agemonths);
cvd_deaths = read.csv(settings$prep.cvd_deaths);
non_cvd_deaths = read.csv(settings$prep.non_cvd_deaths);
hosp_events = read.csv(settings$prep.hosp_events);
assign =  read.csv(settings$prep.assign);
target = readRDS(settings$prep.target)

####################################### CVD subsets #######################################

events = hosp_events
events = subset(events, Incident==1 & Event_Type != "GP" & Event_Type != "Death") # 2265

main_ds = assign %>%
  left_join(agesex, c("id" = "id")) %>%
  left_join(target, c("id" = "Sample_Name")) %>%
  left_join(non_cvd_deaths, c("id" = "id")) %>%
  left_join(cvd_deaths, c("id" = "id"))

ds = main_ds %>% left_join(events, c("id" = "id"))

ds = ds[c("id", "Sample_Sentrix_ID", "sex.x", "age.x", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "Set", "yob", "mob", "dod_ym", 
          "Incident.x", "Event_Type.x", "Event_Date.x", "Incident.y", "Event_Type.y", "Event_Date.y")] 
colnames(ds) = c("id", "Sentrix_ID", "sex", "age", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "set", "yob", "mob", "dod_ym", 
                 "Incident.Death", "Event_Type.Death", "Event_Date.Death", "Incident.Event", "Event_Type.Event", "Event_Date.Event")


####################################### Cox dataset #######################################

cox = make_cox_dataset(ds)
sum(cox$event) # 1831
cox = cox[c("id", "Sentrix_ID", "sex", "age", "assign", "set", "event", "tte")] 
cox = subset(cox, !is.na(tte) & tte>0)
sum(cox$event) # 1825
cox = subset(cox, !is.na(assign)) 
sum(cox$event) # 1624
cox = na.omit(cox)
sum(cox$event) # 1624


####################################### Cox dataset #######################################
boxplot(cox$assign)
sds = 3*sd(cox$assign)
cox_no_nas = subset(cox, !is.na(assign))
outlier_threshold = 1.5*IQR(cox_no_nas$assign) + quantile(cox_no_nas$assign, 0.75)


hist(cox$assign , breaks=40 , col="dodgerblue4" , border=F , main="", xlab="ASSIGN")
dim(subset(cox$assign, assign >= 20))


plot(cox$assign, cox$age)
plot(transform(cox$assign), cox$age)
ggplot(cox, aes(x=age, y=transform(assign))) +
  geom_point() +
  geom_smooth(method='lm', colour = "red") +
  geom_density_2d() +
  xlab("age") +
  ylab("transform(ASSIGN)") +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = 0.5); 

####################################### Basic transformations #######################################

ped = "/Volumes/marioni-lab/Ola/Source/cox/kinship_matrix_using_fixed_2022-01-28-pedigree.rds"
kin_model = readRDS(ped)

cox$assign_transformed = transform(cox$assign)

mod_assign = coxme(Surv(tte, event) ~ scale(age) + factor(sex) + scale(assign_transformed) + (1|id), varlist = kin_model*2, data=cox)
summary(mod_assign)
extract_coxme_table(mod_assign)

between_30_and_70 = cox[cox$age>=30 & cox$age<70, ] # 12518 obs
under_30 = cox[cox$age<30, ]
table(under_30$event)
above_70 = cox[cox$age>=70, ]
table(above_70$event)

mod = coxph(Surv(tte, event) ~ age + sex + scale(assign_transformed), data = cox[cox$age>=30 & cox$age<70, ]) #1352 events
summary(mod)
cox.zph(mod)

mod = coxph(Surv(tte, event) ~ scale(assign_transformed), data = cox[cox$age>=30 & cox$age<70, ]) #1352 events
summary(mod)
cox.zph(mod)

cox$age_cat <- sapply(1:nrow(cox), function(i) {
  if (cox$age[[i]] < 30) {
   return(1)
  } else if (cox$age[[i]]>=30 & cox$age[[i]]<40) {
   return(2)
  } else if (cox$age[[i]]>=40 & cox$age[[i]]<50) {
   return(3)
  } else if (cox$age[[i]]>=50 & cox$age[[i]]<60) {
   return(4)
  } else if (cox$age[[i]]>=60 & cox$age[[i]]<70) {
   return(5)
  } else {
   return(6)
  }
})

df = data.frame("age_category" = numeric(), "HR_assign" = numeric(), "lci" = numeric(), 
                "uci" = numeric(), "p_assign" = numeric(), "n_cases" = numeric(), "n_controls" = numeric())
for (i in 1:6) {
  print(i)
  cat = cox[cox$age_cat==i,]

  mod_assign = coxme(Surv(tte, event) ~ scale(assign_transformed) + (1|id), varlist = kin_model*2, cat)
  coefs = extract_coxme_table(mod_assign)
  mod_asa = coxph(Surv(tte, event) ~ scale(age) + factor(sex) + scale(assign_transformed), data=cat)
  
  tmp_df = data.frame("age_category" = numeric(), "HR_assign" = numeric(), "lci" = numeric(), 
                      "uci" = numeric(), "p_assign" = numeric(), "n_cases" = numeric(), "n_controls" = numeric())
  if(i == 1) {
    tmp_df[1,]$age_category = "<30"
  } else if (i == 6) {
    tmp_df[1,]$age_category = ">=70"
  } else {
    tmp_df[1,]$age_category = paste0((10+i*10), " - ", (20+i*10))
  }
  
  tmp_df[1,]$HR_assign = signif(coefs$hr,3)
  tmp_df[1,]$p_assign = signif(coefs$p,3)
  tmp_df[1,]$lci = signif(coefs$lci,3)
  tmp_df[1,]$uci = signif(coefs$uci,3)
  tmp_df[1,]$n_cases = nrow(subset(cat, event==1))
  tmp_df[1,]$n_controls = nrow(subset(cat, event==0))
  
  if(i == 1) {
    df = tmp_df
  } else {
    df = rbind(df, tmp_df)
  }
}

write.csv(df, '/Volumes/marioni-lab/Ola/Lab/ASSIGN/ASSIGN_HR_split_by_decade_of_age.csv', row.names = F)

ggplot(df, aes(x=HR_assign, y=HR_asa)) +
  geom_point() +
  geom_smooth(method='lm', colour = "red") +
  geom_density_2d() +
  xlab("log(HR_assign)") +
  ylab("log(HR_asa)") +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = 0.5); 
mod = summary(lm(HR.x ~ HR.y, data = merged))
error = calculate_RMSE(mod)

ggplot(df, aes(x=log(HR_assign), y=log(HR_asa))) +
  geom_point() +
  geom_smooth(method='lm', colour = "red") +
  geom_density_2d() +
  xlab("log(HR_assign)") +
  ylab("log(HR_asa)") +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = -1); 
mod = summary(lm(HR.x ~ HR.y, data = merged))
error = calculate_RMSE(mod)

boxplot(cox$assign~cox$event)
sufferers = subset(cox, event == 1)
summary(sufferers$assign)
upper_whisker = as.numeric(quantile(sufferers$assign)[4] + IQR(sufferers$assign) * 1.5) #70
dim(subset(cox, assign>70)) #78

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

cox_old = subset(cox, age>=30 & age<70 & assign < 45) #12665
dim(cox_old)
cox2 = subset(cox, age>=30 & age<70)

boxplot(cox2$assign)
cox2$assign_trimmed = outlierTrim(cox2$assign, 3) 
cox2 = na.omit(cox2) 
dim(cox2)
summary(cox2)
cox2$assign_transformed = transform(cox2$assign)
cox.zph(coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox2))

#outlier threshold by IQR
outlier_threshold_30_70 = 1.5*IQR(cox2$assign) + quantile(cox2$assign, 0.75) #40 
cox_iqr = subset(cox, age>=30 & age<70 & assign <= outlier_threshold_30_70) #12633
cox_iqr$assign_transformed = transform(cox_iqr$assign)
cox.zph(coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox_iqr))

write.csv(cox2[1:8], '/Volumes/marioni-lab/Ola/Lab/Cox/assign_dataset_assumptions_met/age_30-70-3DS-dead-included/age_30-70-3DS-dead-included-12790.csv', row.names = F)             

cox3 = subset(cox, age>=40 & age<70)
cox3$assign_trimmed = outlierTrim(cox3$assign, 3) 
cox3 = na.omit(cox3) 
dim(cox3)
cox3$assign_transformed = transform(cox3$assign)
cox.zph(coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox3))

write.csv(cox3[1:8], '/Volumes/marioni-lab/Ola/Lab/Cox/assign_dataset_assumptions_met/age_40-70-3DS-dead-included/age_40-70-3DS-dead-included-10322.csv', row.names = F)     
          
for(i in 50:60){
  print(i)
  print(cox.zph(coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox[cox$age>=30 & cox$age<70 & cox$assign<i,])))
  print(summary(coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox[cox$age>=30 & cox$age<70 & cox$assign<i,]))$coefficients[1,c(2,5)])
}

cox_filtered = subset(cox, assign<50 & cox$age>=40 & cox$age<70) #n= 12737, number of events= 1270 
mod = coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox_filtered)
summary(mod)
cox.zph(mod)

cox_filtered = subset(cox, assign<45 & cox$age>=30 & cox$age<70) #n= 12245, number of events= 1245 
mod = coxph(Surv(tte, event) ~ scale(assign_transformed), data=cox_filtered)
summary(mod)
cox.zph(mod)

ggcoxfunctional(mod, ylim = c(-0.5,1))
martingale = residuals(mod, type="martingale")

ggcoxdiagnostics(mod, type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(mod, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxzph(cox.zph(mod))

cox$at_risk = 0
cox[which(cox$assign>=20), "at_risk"] = 1

mod1 <- survfit(Surv(tte, event) ~ factor(at_risk), data = cox)
ggsurvplot(mod1,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=FALSE, # show a risk table below the plot
           legend.labs=c("N", "Y"), # change group labels
           legend.title="ASSIGN ≥ 20",  # add legend title
           palette=c("dodgerblue4", "orange"), # change colors of the groups
           ylab="CVD free survival",
           xlab="Time (years)",
           risk.table.height=.2)

