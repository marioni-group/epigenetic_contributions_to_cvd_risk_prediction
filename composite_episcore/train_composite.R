library("optparse")
library("glmnet")
library("survival")
library("jsonlite")

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

# Get arguments
######################################################

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

# folds <- read.delim(opt$folds, header = TRUE, row.names = 2)
df <- read.csv(paste0(settings$out_dir, "cox_covars_episcores_W3_W4.csv"), check.names = FALSE) # 6936 # IQR 6856

# Remove people that have missing or strange time-to-event values (negative)
df = df[!is.na(df$tte) & df$tte>0, ]
df = na.omit(df)

table(df$event)
# 0    1 
# 6275  661  
# IQR
# 0    1 
# 6219  637
# Troponin
# > table(df$event)
# 
# 0    1 
# 6222  658 


if (settings$transform == "rank") {
  df$assign = transform(df$assign)
} else if (settings$transform == "log") {
  df$assign = log(df$assign + 1)
}

x = subset(df, select = -c(Sentrix_ID, id, age, cTnI_corrected, sex, assign, set, event, tte))
x = as.matrix(x)
y = Surv(df$tte,df$event) # Time to event, and wether the event has happened or not
cat("Imported and prepped data!\n")


# Get CV'd lambda and obtain effects of each CpG
######################################################

cat("Fitting elastic net model (Cox PH) for time to CVD...\n")
# Obtain lambda
cv <- cv.glmnet(x, y, family = "cox", type.measure = "C", seed = seed, nfolds = 10) # Harrell's C index to obtain best parameters (kind of like residuals to minimize in OLS)
lambda <- cv$lambda.min
# Obtain coefs
fit <- glmnet(x, y, family = "cox", lambda = lambda, alpha = 0.5)

# Get the good stuff!
######################################################

cat("Now extracting info of interest.\n")
coefs <- coef(fit) # Get betas
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
coefs["Variable"] <- rownames(coefs)
names(coefs)[1] <- "Coefficient"
coefs <- coefs[,c("Variable", "Coefficient")] # 52 coefs # IQR 36 coefs

# Export 
######################################################

cat("Exporting!\n")
filename <- paste0(settings$out_dir, "elnet_coefficients_", settings$run, ".csv")
#filename <- paste0("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_train/w1w3/random/", "elnet_coefficients_random", ".tsv")
write.csv(coefs, file = filename, row.names = FALSE, quote = FALSE)
settings
