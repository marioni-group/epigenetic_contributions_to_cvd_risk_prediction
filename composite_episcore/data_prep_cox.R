library(dplyr)
library(optparse)
library(jsonlite)

#box::use(../../modules/cox[...])
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

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/final/cox_settings.json'

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
assign =  read.csv(settings$prep.assign);
target = readRDS(settings$prep.target)

####################################### CVD subsets #######################################

events = hosp_events
events = subset(events, Incident==1 & Event_Type != "GP" & Event_Type != "Death") # 1297

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

cox = make_cox_dataset(ds) # 18413

# cox = cox[!(cox$dead == 1 & cox$event == 0), ] # remove those who died of non-cardio causes
# cox = cox[c("id", "Sentrix_ID", "sex", "age", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "set", "dead", "event", "tte")] 
cox = cox[c("id", "Sentrix_ID", "sex", "age", "assign", "set", "cTnI_corrected", "event", "tte")] 
cox = subset(cox, !is.na(tte) & tte>0) #18359
cox = subset(cox, !is.na(assign)) # 16354

if(settings$prep.outliers == "SD") {
  cox = subset(cox, cox$age>=30 & cox$age<70) # 12971 items
  cox$assign = outlierTrim(cox$assign, 3) 
  cox = na.omit(cox) #12657
  
  wave1 = subset(cox, set=='wave1')
  dim(wave1)
  table(wave1$event)
  # [1] 3240    9
  # > table(wave1$event)
  # 
  # 0    1 
  # 3240  419 
  
  wave34 = subset(cox, set!='wave1')
  dim(wave34)
  table(wave34$event)
  # > dim(wave34)
  # [1] 8998    9
  # > table(wave34$event)
  # 
  # 0    1 
  # 8143  855 
}
if(settings$prep.outliers == "IQR") {
  outlier_threshold = 1.5*IQR(cox$assign) + quantile(cox$assign, 0.75) # 43
  cox = subset(cox, cox$age>=30 & cox$age<70 & assign <= outlier_threshold) #12633
  
  wave1 = subset(cox, set=='wave1')
  dim(wave1)
  # [1] 3659    9
  > table(wave1$event)
  # 0    1 
  # 3240  419 
  
  wave34 = subset(cox, set!='wave1')
  dim(wave34)
  # [1] 8998    9
  table(wave34$event)
  # 
  # 0    1 
  # 8143  855 
} else {
  # no thresholding
  cox = subset(cox, cox$age>=30 & cox$age<70) #12971
}

write.csv(cox, paste0(settings$out_dir, "cox_covars.csv"), row.names = F)
