#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# written by Elena and Ola

library("optparse")
#library("glmnet")
#library("tidyverse")
library("biglasso")
library("bigmemoryExt")
library("matrixStats")
library("jsonlite")
library("dplyr")

# Note: be careful about excess storage here: /dev/shm!
# Also do lsof /dev/shm for hidden things
# To kill hidden processes: pkill -U ebernab3

# Get arguments
######################################################

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#url = '/Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/settings/sys_bp.json'
url = '/Cluster_Filespace/Marioni_Group/Ola/Code/general/toolbox/troponin_episcores/episcores_using_biglasso/settings/cTnI_corrected_outliers.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)

sink(settings$log_train, split=TRUE)


set.seed(1234) # Set seed to ensure fold variation minimised 
seed <- 1234


# Function for faster row scaling
# Creds: https://www.r-bloggers.com/2016/02/a-faster-scale-function/
######################################################

rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  cm = rowMeans(x, na.rm = TRUE)
  if (scale) {
    csd = rowSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}


# Altered cbindBM for big, big matrices
# OG stuff by cdeterman on GitHub!
######################################################

cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}


# x needs to be a list of big matrices, and they have to have the same rownames
cbindBM_list <- function(x, binding="right", 
                         z=NULL, type=NULL, separated=NULL,
                         backingfile=NULL, backingpath=NULL,
                         descriptorfile=NULL, binarydescriptor=FALSE,
                         shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
    cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
    cols_list <- append(cols_list, list(cols))
    total_cols <- total_cols + ncol(x[[i]])
  }    
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
      z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
      z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)
    
    if (erase == TRUE) {
      cat("\nErasing chunk and liberating memory...\n\n")
      x[[i]] <- "Replacement"
      gc()
    }
  }
  return(z)
}

lasso <- settings$lasso

pheno = read.csv(settings$o_pheno_train)
pheno = subset(pheno, !duplicated(pheno$Sample_Sentrix_ID))
rownames(pheno) = pheno$Sample_Sentrix_ID
cat("\nImported pheno + fold data!\n")
cat(paste0("Remaining individuals: ", length(rownames(pheno)),"\n"))

meth <- readRDS(settings$o_train_df)
meth = subset(meth, !duplicated(rownames(meth)))
meth <- meth[which(rownames(meth) %in% rownames(pheno)),]
pheno <- pheno[which(rownames(pheno) %in% rownames(meth)),] 
cat("\nImported methylation data!\n")

# Import methylation plus pheno data
######################################################


cat("\nData prep...\n")


# Check order and divide if sex-stratifying
######################################################

# Check rownames are the same and in the same order

if (identical(rownames(meth), rownames(pheno))) {
  cat("\nRownames match.\n")
} else {
  meth <- meth[match(rownames(pheno), rownames(meth)),]
  cat("\nRownames have been matched.\n")
}

cat("\nRAM clean up...\n\n")
gc()



# Create big matrix object for biglasso
######################################################

cat("\nMaking big matrix objects for big lasso...\n")
div <- 5 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()
  
for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
x <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)


cat("\nRAM clean up...\n\n")
gc()


# Stuff for biglasso
######################################################

if (length(rownames(pheno)) < 10000) { #w1w3
  nfolds <- 10
} else if ( (length(rownames(pheno)) > 10000) & (length(rownames(pheno)) < 15000)) { #w1w3 + LBC + GEO
  nfolds <- 15
} else if ( (length(rownames(pheno)) > 15000) & (length(rownames(pheno)) < 20000)) { #w1w3w4
  nfolds <- 20
} else if (length(rownames(pheno)) > 20000) { #w1w3w4 + LBC + GEO
  nfolds <- 25
}

cat(paste0("\nTotal number of folds: ", nfolds,"\n"))
folds <- pheno$Fold
y <- as.numeric(pheno[[settings$feature]])

cat(paste0("x dimensions:\n"))
cat(dim(x))
cat("\n")
cat(paste0("y dimensions:\n"))
cat(length(y))

# Lasso/elnet parameters
penalty <- "enet"
alpha <- 0.5
if (lasso == "T") {
  penalty <- "lasso"
  alpha <- 1
}

cat("\n\nRAM clean up...\n\n")
gc()


# Get CV'd lambda and obtain effects of each CpG
######################################################

cat("\nFitting elastic net model for age...\n")
cat("\nMode selected: randomized folds.\n")
cvfit <- cv.biglasso(x, y, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, nfolds = nfolds, penalty = penalty)

lambda <- cvfit$lambda.min # Get lambda which gives minimum MSE
# Obtain coeffs
fit <- biglasso(x, y, family = "gaussian", alpha = alpha, ncores = 8, lambda = lambda, penalty = penalty)
cat("\nRAM clean up...\n\n")
gc()


# Get the good stuff!
######################################################

cat("\nNow extracting info of interest.\n")
coefs <- coef(fit) # Get betas
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient"
coefs["CpG_loc"] <- noquote(sub('V', '', rownames(coefs)))
coefs[1, "CpG_loc"] <- NA
coefs[["CpG_loc"]] <- as.numeric(coefs[["CpG_loc"]])
coefs[2:nrow(coefs), "Variable"] <- names[coefs[2:nrow(coefs), "CpG_loc"]]
coefs[1, "Variable"] <- "Intercept"
coefs <- coefs[c("Variable", "Coefficient")]



# Export coeffs
######################################################

cat("\nExporting!\n\n")
filename <- paste0(settings$o_dir, "elnet_coefficients_", settings$feature, ".csv")
#filename <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/elnet_train/w1w3w4/random_noscalesample_noscalecpg_subset20K/elnet_coefficients_random_noscalesample_noscalecpg_subset20K_squaredsubset900.tsv"
write.csv(coefs, file = filename, row.names = FALSE)

sink()

