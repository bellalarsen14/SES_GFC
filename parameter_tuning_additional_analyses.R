# load libraries
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

root <- "/parentfoldername/"

# define the behavioral (SES) variables for analysis
behavvar_list <- data.frame(behavvar = c("ADI311","SESchildhd",
                                         "PH45_AreaDeptot","SESall45",
                                         "neighdep2645_factor","ses_composite"))

# define the reliability threshold for edges at 0.75
ICCthr <- 0.75

# define the brain-based variable for analysis: general functional connectivity (GFC)
brainvarlists <- list(list("GFC")) 

# define the number of iterations for tuning
n_iter <- 100

# set up loop configurations
  # behavvar_idx / covariate_idx index into behavvar_list (1-6)
  # filter_var1 / filter_var2 are the two variables used to restrict to
  # complete cases before splitting (within SES level, consistent from childhood to adulthood)
run_configs <- data.frame(
  name          = c("adult_neigh", "adult_ses"),
  behavvar_idx  = c(5, 6),
  covariate_idx = c(1, 2),
  filter_var1   = c("neighdep2645_factor", "ses_composite"),
  filter_var2   = c("ADI311", "SESchildhd"),
  workdir1_name = c("adult_neigh_covar_child", "adult_ses_covar_child"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Read data 
# ---------------------------------------------------------------------------

## load FC edges
### loads variable "FC_ALL", columns are "snum" and "edge#"
load(paste0(root, 'DBIS_GFC_N769_incSubcortex.Rdata'))
GFC <- FC_ALL
GFC$id <- as.numeric(sub("sub-","",GFC$id)); names(GFC)[1] <- "snum"

## ROIs
ROIs_GFC <- names(GFC)[grepl("edge", names(GFC))]

## load ICC values (1 value per edge, vector of 78210 ICC values)
ICCs_GFC <- read.csv(paste0(root,'DBIS_GFC_N769_incSubcortex_ICCs.csv'))$ICC

## load behavioral data file, a dataframe with subject number, one row per participant, and SES data
behavdata <- read.csv(file = paste0(root,"person_level_df_all_ses_comp.csv")) %>%
  select(-X)

## ensure that sex is a factor with levels 1,2 where 1 is the reference level
behavdata$sex <- as.factor(behavdata$sex)

## load in motion (framewise displacement, or FD) covariate, dataframe with subject number, one row per participant, and FD value
motion <- read.csv(paste0(root,'Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_motion.csv'))
motion$snum <- as.numeric(sub("sub-","",motion$id))

## join behavioral data with imaging data to pare down to subjects with both data
behav_merged_full <- join_all(list(
  behavdata,
  motion[,c("snum","AverageFD")]
), by="snum", type="full")

# make GFC dataset 
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

  workdir1 <- paste0(root, "Larsen/Dunedin/Bella_Prediction_Outputs/Tuning_loop/", cfg$workdir1_name, "/enet") # need to create empty storage folders before running
  
  # pair-specific complete-case filter
  behav_merged <- behav_merged_full %>%
    filter(!is.na(.data[[cfg$filter_var1]]) & !is.na(.data[[cfg$filter_var2]]))

  # create dataframe for base model, with sex and motion as covariates
  data <- dplyr::left_join(braindat_scaled[,c("snum",ROIs_full)], behav_merged[, c("snum","sex","AverageFD",behavvar)], by="snum")
  data <- data[complete.cases(data),]

  # create dataframe for covariate model, with sex, motion, and alternate SES level as covariates
  data_covar <- dplyr::left_join(braindat_scaled[,c("snum",ROIs_full)], behav_merged[, c("snum","sex","AverageFD",behavvar,covariate_1)], by="snum")
  data_covar <- data_covar[complete.cases(data_covar),]

  runname <- 'test'

  # -------------------------------------------------------------------------
  # Pre-generate all n_iter train/test splits up front, before any model
  # fitting. Ensures consistency across runs between the models with and 
  # without covariates.
  # -------------------------------------------------------------------------
  set.seed(12345)
  splits_cur <- lapply(1:n_iter, function(i) {
    createDataPartition(data[, paste(behavvar)], p = 0.9, list = FALSE)
  })

  # save splits so other scripts (e.g. no-covariate version) can load
  # the exact same object instead of regenerating them.
  assign(paste0("splits_", cfg$name), splits_cur)
  save(list = paste0("splits_", cfg$name),
       file = paste0(root, "Tuning/generated_splits_loop/adult_covar_chld/splits_", cfg$name, ".Rdata"))

  for (iter in 1:n_iter){
    perf <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), RMSE=numeric(), Rsquare=numeric(), r=numeric(), MAE=numeric())

    # use pre-generated test/train split for this iteration
    # note: a test/train split is created here though only training data is
    # used for tuning; this choice was made to test robustness of parameter
    # choice without creating leakage. This retains blindness to the test
    # data. Lambda is tuned within train.data only.
    training.samples <- splits_cur[[iter]]
    train.data <- data[training.samples, ]
    test.data  <- data[-training.samples, ]

    train_snums <- train.data$snum
    test_snums <- test.data$snum
    
     # add covariate and run
    perf_covar <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), RMSE=numeric(), Rsquare=numeric(), r=numeric(), MAE=numeric())
    
    # use the same test/train datasets as previously defined w/o covariates but rename to avoid writing over
    train.data.covar <- data_covar[training.samples, ]
    test.data.covar  <- data_covar[-training.samples, ]
    
    # regress sex and motion, and covariate_1, from training set
    lm_covar <- lm(train.data.covar[,paste(behavvar)] ~ train.data.covar$sex + train.data.covar$AverageFD + train.data.covar[,paste(covariate_1)])
    train.data.covar$behav_resids <- scale(lm_covar$residuals) #### sure we want to scale this way?
    # adjust using same parameters in test set
    coefs_covar <- lm_covar$coefficients # 1=intercept, 2=sex, 3=averageFD, 4=covariate
    test_fitted_covars <- coefs_covar[1] + coefs_covar[2]*(as.numeric(test.data.covar$sex)-1) +
      coefs_covar[3]*test.data.covar$AverageFD + coefs_covar[4]*test.data.covar[,paste(covariate_1)]
    test.data.covar$behav_adj <- scale( test.data.covar[,paste(behavvar)] - test_fitted_covars )
    
    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    train.data.covar <- train.data.covar[, c("behav_resids", ROIs)]
    test.data.covar <- test.data.covar[, c("behav_adj", ROIs)]
    
    # Setup a grid range of lambda values:
    lambdas <- 10^seq(-2, 2, length = 25)
    
    ## predict with ridge regression
    ridge_covar <- train(
      as.formula(paste("behav_resids", "~ .")), data = train.data.covar, method = "glmnet",
      #trControl = ctrl_frozen,
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambdas)
    )
    
    enet_covar <- ridge_covar
    
    ## save out everything
    outname <- paste0(gsub(" ", "_", gsub(":","_",date())), "_", round(runif(1,100,999),0))
    save(enet_covar, file=paste0(workdir1,"/tuning_",outname,".Rdata"))

  }
}
