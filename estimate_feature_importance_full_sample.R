# adapted from predictions_cluster scripts in gradientreliabilitypaper folder for adding in cerebellum
library(haven)
library(readr)
library(caret)
library("matrixStats")
library(tidyr)
library(plyr)
library(foreign)
library(dplyr)
library(glmnet)
library(glmnetUtils) # for more memory-efficient cv ridge regression. must be loaded after glmnet!
suppressMessages(require(optparse))

################### OPTIONS #############################

root <- "parentfoldername"
behavvar_list <- read.csv(paste0(root,'behavvar_list.csv'))

ICCthr <- 0.75
n_iter <- 1 # run only once in the full sample

# each brainvarlist is a list of data frame names. They should have corresponding variables for ROIs within the DF, named ROIs_<DF>
# brainvarlists <- list(list("asegALL", "CT", "SA", "GMV")) 
# brainvarlists <- list(list("asegALL", "gCT", "gSA", "gGMV")) 
# brainvarlists <- list(list("asegALL", "gCT", "gSA", "gGMV", "FA","MD","AD","RD","WMH","GFC")) 
brainvarlists <- list(list("GFC")) 

################### CONFIG: one row per variable/pair ####################
# behavvar_idx / covariate_idx index into behavvar_list (1-6)
# filter_var1 / filter_var2 are the two variables used to restrict to
#   complete cases before splitting (identical for both members of a pair)
# split_obj is the name of the object loaded from each pair's saved splits file
run_configs <- data.frame(
  name           = c("chldhd_neigh", "chldhd_ses", "adult_neigh", "adult_ses", "age45_neigh", "age45_ses"),
  behavvar_idx   = c(1, 2, 5, 6, 3, 4),
  covariate_idx  = c(2, 1, 6, 5, 4, 3),
  filter_var1    = c("SESchildhd", "SESchildhd", "neighdep2645_factor", "neighdep2645_factor", "PH45_AreaDeptot", "PH45_AreaDeptot"),
  filter_var2    = c("ADI311", "ADI311", "ses_composite", "ses_composite", "SESall45", "SESall45"),
  folder_name    = c("chldhd_neigh_full", "chldhd_ses_full", "adult_neigh_full", "adult_ses_full", "age45_neigh_full", "age45_ses_full"),
  stringsAsFactors = FALSE
)

################### END OPTIONS #############################

# ---------------------------------------------------------------------------
# Read data 
# ---------------------------------------------------------------------------

## load FC edges - takes about a minute! (using the function to load individually took over an hour, so finally just saved out the results!)
### loads variable "FC_ALL", columns are "snum" and "edge#"
load(paste0(root, 'DBIS_GFC_N769_incSubcortex.Rdata'))
GFC <- FC_ALL
GFC$id <- as.numeric(sub("sub-","",GFC$id)); names(GFC)[1] <- "snum"

## ROIs
ROIs_GFC <- names(GFC)[grepl("edge", names(GFC))]

## define the reliability threshold for edges at 0.75
ICCs_GFC <- read.csv(paste0(root,'DBIS_GFC_N769_incSubcortex_ICCs.csv'))$ICC

# load behavioral data file, a dataframe with subject number, one row per participant, and SES data
behavdata <- read.csv(file = paste0(root,"person_level_df_all_ses_comp.csv")) %>%
  select(-X)

# make sex a factor
behavdata$sex <- as.factor(behavdata$sex)

# covariates
motion <- read.csv(paste0(root,'Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_motion.csv'))
motion$snum <- as.numeric(sub("sub-","",motion$id))

# NOTE: unfiltered join -- the pair-specific complete-case filter
# (filter_var1 / filter_var2) is applied per-config inside the loop.
behav_merged_full <- join_all(list(
  behavdata,
  motion[,c("snum","AverageFD")]
), by="snum", type="full")

# make GFC dataset (ROIs_full / braindat_scaled do not depend on config,
# since brainvarlists and ICCthr are fixed across all 6 runs)
for ( brainvarlist in brainvarlists ) {
  braindat <- data.frame(snum=GFC$snum)
  ROIs_full <- c()
  for (i in 1:length(brainvarlist)) {
    brainvarcur <- brainvarlist[[i]]
    ROIs_cur <- get(paste0("ROIs_", brainvarcur))
    ICCs_cur <- get(paste0("ICCs_", brainvarcur))
    # threshold ROIs by ICC
    ROIs_full <- c(ROIs_full, ROIs_cur[ICCs_cur > ICCthr])
    braindat <- merge(braindat, get(brainvarcur), by="snum")
    #Fisher's Z adjust the ROIs
    braindat_scaled <- braindat
    braindat_scaled[, -1] <- atanh(braindat[, -1])
    if (i>1) { brainvar <- paste(brainvar, brainvarcur, sep=".") } else { brainvar <- brainvarcur }
  }
}

# ---------------------------------------------------------------------------
# Main loop: one iteration per config row
# ---------------------------------------------------------------------------
for (cfg_row in 1:nrow(run_configs)) {
  
  cfg <- run_configs[cfg_row, ]
  
  behavvar <- behavvar_list[cfg$behavvar_idx, ]
  covariate_1 <- behavvar_list[cfg$covariate_idx, ]
  
  workdir1 <- paste0(root, "Updated_Runs_CV/", cfg$folder_name, "/haufe_coef_nocovar")
  workdir4 <- paste0(root, "Updated_Runs_CV/", cfg$folder_name, "/haufe_coef_covar")
  #create working directory
  for (wd in c(workdir1, workdir4)) dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  
  # pair-specific complete-case filter
  behav_merged <- behav_merged_full %>%
    filter(!is.na(.data[[cfg$filter_var1]]) & !is.na(.data[[cfg$filter_var2]]))
  
  data <- dplyr::left_join(braindat_scaled[,c("snum",ROIs_full)], behav_merged[, c("snum","sex","AverageFD",behavvar)], by="snum")
  data <- data[complete.cases(data),]
  
  data_covar <- dplyr::left_join(braindat_scaled[,c("snum",ROIs_full)], behav_merged[, c("snum","sex","AverageFD",behavvar,covariate_1)], by="snum")
  data_covar <- data_covar[complete.cases(data_covar),]
  
  runname <- 'test'
  
  
  for (iter in 1:n_iter){

    # Fix the CV folds used to tune lambda so the Haufe coefficients are reproducible.
    set.seed(98765 + cfg_row)
    
    perf <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), RMSE=numeric(), Rsquare=numeric(), r=numeric(), MAE=numeric())
    
    # this time, we are running the model only once, on the full data sample
    
    # regress sex and motion from full data set
    lm <- lm(data[,paste(behavvar)] ~ data$sex + data$AverageFD)
    data$behav_resids <- scale(lm$residuals) 

    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    data <- data[, c("behav_resids", ROIs)]
    
    # Setup a grid range of lambda values:
    lambdas <- 10^seq(-2, 2, length = 25)
    
    ## train with ridge regression
    ridge <- train(
      as.formula(paste("behav_resids", "~ .")), data = data, method = "glmnet",
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambdas)
    )
    
    ## predict once per variable in FULL data
    predictions_ridge_train <- ridge %>% predict(data)
    
    ## haufe transform for coefficients per Tian and Zalesky NI 2021
    coefs_haufe <- c()
    N <- nrow(data)
    for (r in ROIs){ # loop through all edges
      r_std <- scale(data[,paste(r)])
      coefs_haufe <- c( coefs_haufe, sum(r_std * predictions_ridge_train) / N )
    }
    
    ## save out everything
    outname <- paste0(gsub(" ", "_", gsub(":","_",date())), "_", round(runif(1,100,999),0))
    save(coefs_haufe, file=paste0(workdir1,"/predictions_",outname,".Rdata"))
    
    # add covariate and run again
    perf_covar <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), RMSE=numeric(), Rsquare=numeric(), r=numeric(), MAE=numeric())
    
    # regress sex and motion, and covariate_1, from training set
    lm_covar <- lm(data_covar[,paste(behavvar)] ~ data_covar$sex + data_covar$AverageFD + data_covar[,paste(covariate_1)])
    data_covar$behav_resids <- scale(lm_covar$residuals) 

    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    data_covar <- data_covar[, c("behav_resids", ROIs)]
    
    ## train with ridge regression
    ridge_covar <- train(
      as.formula(paste("behav_resids", "~ .")), data = data_covar, method = "glmnet",
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambdas)
    )
    
    #   # haufe transform for coefficients per Tian and Zalesky NI 2021
    predictions_ridge_train_covar <- ridge_covar %>% predict(data_covar)
    coefs_haufe_covar <- c()
    N_covar <- nrow(data_covar)
    for (r in ROIs){ # loop through all edges
      r_std <- scale(data_covar[,paste(r)])
      coefs_haufe_covar <- c( coefs_haufe_covar, sum(r_std * predictions_ridge_train_covar) / N_covar )
    }
    
    ## save out everything
    outname <- paste0(gsub(" ", "_", gsub(":","_",date())), "_", round(runif(1,100,999),0))
    save(coefs_haufe_covar, file=paste0(workdir4,"/predictions_",outname,".Rdata"))
  }
}

