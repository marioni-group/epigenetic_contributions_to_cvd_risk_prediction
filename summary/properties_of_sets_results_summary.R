library(dplyr)
library(optparse)
library(jsonlite)
library(survival)
library(survminer)

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

make_cox_dataset = function(ds) {
  
  ds[which(is.na(ds$mob)), "mob"] = 1
  
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

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
path = ''
url = paste0(path, 'Ola/Lab/EpiScores/Cox_episcores_composite/runs/final/cox_settings.json')

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)
settings

set.seed(42) # Set seed to ensure fold variation minimised 
seed <- 42


####################################### Input data #######################################

agesex = read.csv(settings$prep.agesex);
agemonths = read.csv(settings$prep.agemonths);
cvd_deaths = read.csv(settings$prep.cvd_deaths);
non_cvd_deaths = read.csv(settings$prep.non_cvd_deaths);
hosp_events = read.csv(settings$prep.hosp_events);
assign =  read.csv(paste0(path, 'Ola/Lab/ASSIGN/GS_assign_data_for_summary.csv'))
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
          "Incident.x", "Event_Type.x", "Event_Date.x", "Incident.y", "Event_Type.y", "Event_Date.y", "simd", "history_cvd", "diabetic", "ra",
          "cigs_day_raw", "sys_bp", "total_cholesterol", "HDL_cholesterol")] 
colnames(ds) = c("id", "Sentrix_ID", "sex", "age", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "set", "yob", "mob", "dod_ym", 
                 "Incident.Death", "Event_Type.Death", "Event_Date.Death", "Incident.Event", "Event_Type.Event", "Event_Date.Event", "simd", "history_cvd", "diabetic", "ra",
                 "cigs_day_raw", "sys_bp", "total_cholesterol", "HDL_cholesterol")


####################################### Cox dataset #######################################

cox = make_cox_dataset(ds)
cox$non_smoker = ifelse(cox$cigs_day_raw==0, 1, 0)
sum(cox$event) # 1831

cox = subset(cox, !is.na(tte) & tte>0)
sum(cox$event) # 1825
cox = subset(cox, !is.na(assign)) # 15664
sum(cox$event) # 1624, 16354
cox = subset(cox, cox$age>=30 & cox$age<70) # 12971 items
cox$assign = outlierTrim(cox$assign, 3) 
cox = subset(cox, !is.na(assign))
table(cox$event) # 1288

write.csv(cox, paste0(path, 'Ola/Lab/Cox/basic_dataset_no_assumptions/cox_covars.csv'), row.names = F) # 15664 items

####################################### By wave #######################################

only_assign_vars = cox[c("id", "Sentrix_ID", "sex", "age", "simd", "history_cvd", "diabetic", 
                         "cigs_day_raw", "non_smoker", "sys_bp", "ra", "total_cholesterol","HDL_cholesterol",
                         "assign", "set", "event", "tte", "Troponin_I", "cTnI_corrected")] 
only_assign_vars = na.omit(only_assign_vars)
table(only_assign_vars$event) # 1288

w1 = subset(only_assign_vars, set=="wave1")
w3_w4 = subset(only_assign_vars, set!="wave1") #9079
w3 = subset(only_assign_vars, set=="wave3")
w4 = subset(only_assign_vars, set=="wave4")

unrelated = read.csv(settings$prep.unique_W3_W4)
w3_w4 = w3_w4[which(w3_w4$id %in% unrelated$Sample_Name),] #6936

table(is.na(w1))
table(is.na(w3_w4))
dim(w1)
dim(w3_w4)

table(w1$event)
table(w3_w4$event)

threshold_events = function(x, threshold = 10) {
  x$event <- sapply(1:nrow(x), function(i) {
    if (x$tte[[i]] > threshold) {
      0
    } else {
      x$event[[i]]
    }
  })
  return(x)
}

w1 = threshold_events(w1, 10)
table(w1$event)

full_cases = subset(only_assign_vars, event == 1)
full_controls = subset(only_assign_vars, event == 0)
training_cases = subset(w3_w4, event == 1)
training_controls = subset(w3_w4, event == 0)
test_cases = subset(w1, event == 1)
test_controls = subset(w1, event == 0)

write.csv(cox, paste0(path, 'Ola/Lab/Cox/basic_dataset_no_assumptions/cox_covars.csv'), row.names = F) # 15664 items

####################################### Subsets #############################################

continous_normal_stats = function(variable, set) {
  mean = mean(set[,variable])
  sd = sd(set[,variable])
  ret = c(mean, sd)
  return(ret)
}

continous_non_normal_stats = function(variable, set) {
  median = median(set[,variable])
  q1 = quantile(set[,variable], 0.25)
  q3 = quantile(set[,variable], 0.75)
  ret = c(median, q1, q3)
  return(ret)
}

categorical_stats = function(variable, response, set) {
  condition = paste0(variable, "==", response)
  tmp_set = subset(set, eval(parse(text=condition)))
  n_set = nrow(set)
  n_tmp_set = nrow(tmp_set)
  prop = n_tmp_set * 100/n_set
  ret = c(n_tmp_set, prop)
  return(ret)
}

gather_set_statistics = function(set, set_name) {
  tte = continous_non_normal_stats("tte", set)
  age = continous_non_normal_stats("age", set)
  sex = categorical_stats("sex", 1, set)
  simd = continous_non_normal_stats("simd", set)
  history = categorical_stats("history_cvd", 1, set)
  diabetes = categorical_stats("diabetic", 1, set)
  ra = categorical_stats("ra", 1, set)
  non_smoker = categorical_stats("non_smoker", 1, set)
  sys_bp = continous_normal_stats("sys_bp", set)
  total_cholesterol = continous_normal_stats("total_cholesterol", set)
  hdl = continous_non_normal_stats("HDL_cholesterol", set)
  assign_score = continous_non_normal_stats("assign", set)
  
  formatted_list = data.frame(
    "set" = set_name,
    "n" = nrow(set),
    "tte" = sprintf("%.1f\r\n[%.1f, %.1f]", tte[1], tte[2], tte[3]),
    "age" = sprintf("%.1f\r\n[%.1f, %.1f]", age[1], age[2], age[3]),
    "sex"= sprintf("%.f\r\n(%.1f%%)", sex[1], sex[2]),
    "simd" = sprintf("%.1f\r\n[%.1f, %.1f]", simd[1], simd[2], simd[3]),
    "history_cvd" = sprintf("%.f\r\n(%.1f%%)", history[1], history[2]),
    "diabetes" = sprintf("%.f\r\n(%.1f%%)", diabetes[1], diabetes[2]),
    "ra" = sprintf("%.f\r\n(%.1f%%)", ra[1], ra[2]),
    "non_smoker" = sprintf("%.f\r\n(%.1f%%)", non_smoker[1], non_smoker[2]),
    "sys_bp" = sprintf("%.1f\r\n(%.1f)", sys_bp[1], sys_bp[2]),
    "total_cholesterol" = sprintf("%.1f\r\n(%.1f)", total_cholesterol[1], total_cholesterol[2]),
    "HDL_cholesterol" = sprintf("%.1f\r\n[%.1f, %.1f]", hdl[1], hdl[2], hdl[3]),
    "ASSIGN" = sprintf("%.f\r\n[%.f, %.f]", assign_score[1], assign_score[2], assign_score[3])
  )
  
  return(formatted_list)
}

sets = list(training_cases = training_cases, 
            training_controls = training_controls, 
            test_cases = test_cases, 
            test_controls = test_controls, 
            full_cases = full_cases, 
            full_controls = full_controls)

statistics <- lapply(seq_along(sets), function(i) {gather_set_statistics(sets[[i]], names(sets)[[i]])})
statisticsAssign = t(bind_rows(statistics))

write.csv(statisticsAssign, paste0(path, 'Ola/Lab/ASSIGN/summaryTable.csv'), row.names = T)
