library(dplyr)
library(parallel)
library(optparse)
library(jsonlite)
box::use(../../modules/transformations[...])

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/episcores_only_valid_assumptions_10_years/composite_cox_settings.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)
settings

set.seed(42) # Set seed to ensure fold variation minimised 
seed <- 42
cox = read.csv(paste0(settings$out_dir, "cox_covars.csv"))
min_set = subset(cox, !is.na(tte) & tte>0)
min_set = subset(min_set, !is.na(assign))

episcores_w1_w3 = read.csv(settings$prep.episcores_W1_W3, check.names = FALSE)
episcores_w4 = read.csv(settings$prep.episcores_W4, check.names = FALSE)
names(episcores_w4)[1] = 'ID'
identical(colnames(episcores_w1_w3), colnames(episcores_w4))
episcores = union(episcores_w1_w3, episcores_w4)
merged = merge(min_set, episcores, by.x="Sentrix_ID", by.y="ID")
# correct assumptions
merged = subset(merged, tte<10)

#join with cox
W1 = subset(merged, set == "wave1")
W3_W4 = subset(merged, set == "wave4" | set == "wave3")

#keep only unrelated
unrelated = read.csv(settings$prep.unique_W3_W4)
W3_W4 = W3_W4[which(W3_W4$id %in% unrelated$Sample_Name),] #8791

start = length(min_set) + 1
end = ncol(merged)

head(W3_W4[,1:10])
head(W3_W4[,100:ncol(W3_W4)])
W3_W4[start:end] = mclapply(W3_W4[start:end], transform)
head(W3_W4[,100:ncol(W3_W4)])
head(W3_W4[,1:10])

head(W1[,1:10])
head(W1[,100:ncol(W1)])
W1[start:end] = mclapply(W1[start:end], transform)
head(W1[,100:ncol(W1)])
head(W1[,1:10])

#check col_length
table(is.na(W1))
table(is.na(W3_W4))

write.csv(W1, paste0(settings$out_dir, "cox_covars_episcores_W1.csv"), row.names = F)
write.csv(W3_W4, paste0(settings$out_dir, "cox_covars_episcores_W3_W4.csv"), row.names = F)
