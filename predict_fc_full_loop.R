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

root <- "/parentfoldername/"

# define the behavioral (SES) variables for analysis
behavvar_list <- data.frame(behavvar = c("ADI311","SESchildhd",
                                         "PH45_AreaDeptot","SESall45",
                                         "neighdep2645_factor","ses_composite"))

# load in hyperparameters selected from tuning 
lambda_selection <- read.csv(paste0(root,'lambda_selection.csv'))

# define the reliability threshold for edges at 0.75
ICCthr <- 0.75

# define the brain-based variable for analysis: general functional connectivity (GFC)
brainvarlists <- list(list("GFC")) 

# define the number of iterations for tuning
n_iter <- 100

# set up loop configurations
  # behavvar_idx / covariate_idx index into behavvar_list (1-6)
  # lambda_row / lambda_covar_row index into lambda_selection (row, then column 2)
  # filter_var1 / filter_var2 are the two variables used to restrict to
  #   complete cases before splitting (identical for both members of a pair)
  # split_obj is the name of the object loaded from each pair's saved splits file (from tuning)
run_configs <- data.frame(
  name           = c("chldhd_neigh", "chldhd_ses", "adult_neigh", "adult_ses", "age45_neigh", "age45_ses"),
  behavvar_idx   = c(1, 2, 5, 6, 3, 4),
  covariate_idx  = c(2, 1, 6, 5, 4, 3),
  lambda_row     = c(10, 12, 2, 4, 6, 8),
  lambda_covar_row = c(9, 11, 1, 3, 5, 7),
  filter_var1    = c("SESchildhd", "SESchildhd", "neighdep2645_factor", "neighdep2645_factor", "PH45_AreaDeptot", "PH45_AreaDeptot"),
  filter_var2    = c("ADI311", "ADI311", "ses_composite", "ses_composite", "SESall45", "SESall45"),
  folder_name    = c("chldhd_neigh_full", "chldhd_ses_full", "adult_neigh_full", "adult_ses_full", "age45_neigh_full", "age45_ses_full"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------

## load FC edges
### loads variable "FC_ALL", columns are "snum" and "edge#"
load(paste0(root, 'Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_incSubcortex.Rdata'))
GFC <- FC_ALL
GFC$id <- as.numeric(sub("sub-","",GFC$id)); names(GFC)[1] <- "snum"

## ROIs
ROIs_GFC <- names(GFC)[grepl("edge", names(GFC))]

## load ICC values (1 value per edge, vector of 78210 ICC values)
ICCs_GFC <- read.csv(paste0(root,'Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_incSubcortex_ICCs.csv'))$ICC

## load behavioral data file, a dataframe with subject number, one row per participant, and SES data
behavdata <- read.csv(file = paste0(root,"Larsen/Dunedin/Bella_Code/person_level_df_all_ses_comp.csv")) %>%
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

  workdir1 <- paste0(root, "Larsen/Dunedin/Bella_Prediction_Outputs/New_Runs_Lambda_Tuned_LOOP/", cfg$folder_name, "/pred_nocovar/") # need to create empty storage folders before running
  workdir2 <- paste0(root, "Larsen/Dunedin/Bella_Prediction_Outputs/New_Runs_Lambda_Tuned_LOOP/", cfg$folder_name, "/perf_nocovar/") # need to create empty storage folders before running
  workdir3 <- paste0(root, "Larsen/Dunedin/Bella_Prediction_Outputs/New_Runs_Lambda_Tuned_LOOP/", cfg$folder_name, "/perf_covar_ses/") # need to create empty storage folders before running
  workdir4 <- paste0(root, "Larsen/Dunedin/Bella_Prediction_Outputs/New_Runs_Lambda_Tuned_LOOP/", cfg$folder_name, "/pred_covar_ses/") # need to create empty storage folders before running

  lambda_set <- lambda_selection[cfg$lambda_row, 2]         # from pre-run tuning
  lambda_set_covar <- lambda_selection[cfg$lambda_covar_row, 2]

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
  # Load the matching pre-generated splits for this config (produced by the
  # tuning script, so the same 100 splits are reused here).
  # -------------------------------------------------------------------------
  load(paste0(root, "Larsen/Dunedin/Code\ from\ Ethan/Bella\ edits/Tuning/generated_splits_loop/splits_", cfg$name, ".Rdata"))  # loads object: splits_<name>
  splits_cur <- get(paste0("splits_", cfg$name))

  for (iter in 1:n_iter){
    perf <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), RMSE=numeric(), Rsquare=numeric(), r=numeric(), MAE=numeric())

    # use pre-generated test/train split for this iteration
    training.samples <- splits_cur[[iter]]
    train.data <- data[training.samples, ]
    test.data  <- data[-training.samples, ]

    train_snums <- train.data$snum
    test_snums <- test.data$snum

    # regress sex and motion from training set
    lm <- lm(train.data[,paste(behavvar)] ~ train.data$sex + train.data$AverageFD)
    train.data$behav_resids <- scale(lm$residuals) #### sure we want to scale this way?
    # adjust using same parameters in test set
    coefs <- lm$coefficients # 1=intercept, 2=sex, 3=averageFD
    test_fitted_covars <- coefs[1] + coefs[2]*(as.numeric(test.data$sex)-1) +
      coefs[3]*test.data$AverageFD
    test.data$behav_adj <- scale( test.data[,paste(behavvar)] - test_fitted_covars )

    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    train.data <- train.data[, c("behav_resids", ROIs)]
    test.data <- test.data[, c("behav_adj", ROIs)]

    ## train with ridge regression
    ridge <- train(
      as.formula(paste("behav_resids", "~ .")), data = train.data, method = "glmnet",
      #trControl = ctrl_frozen,
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambda_set)
    )
    
    ## predict in test data
    predictions_ridge <- ridge %>% predict(test.data)
    perf <- data.frame( method = "ridge", brainvar = 'GFC', behavvar=behavvar, nROIs = length(ROIs), iteration = iter, N=nrow(data),
                        RMSE = RMSE(predictions_ridge, test.data$behav_adj),
                        Rsquare = R2(predictions_ridge, test.data$behav_adj),
                        MAE = MAE(predictions_ridge, test.data$behav_adj),
                        r = cor(predictions_ridge, test.data$behav_adj ) )

    ## haufe transform for coefficients per Tian and Zalesky NI 2021
    predictions_ridge_train <- ridge %>% predict(train.data)
    coefs_haufe <- c()
    N <- nrow(train.data)
    for (r in ROIs){ # loop through all edges
      r_std <- scale(train.data[,paste(r)])
      coefs_haufe <- c( coefs_haufe, sum(r_std * predictions_ridge_train) / N )
    }

    ## save out everything
    outname <- paste0(gsub(" ", "_", gsub(":","_",date())), "_", round(runif(1,100,999),0))
    df <- data.frame(snum=test_snums, prediction_ridge=predictions_ridge)
    save(df, coefs_haufe, ROIs_full, file=paste0(workdir1,"/predictions_",outname,".Rdata"))
    save(perf, file=paste0(workdir2,"/performance_",outname,".Rdata"))

    # add covariate and run again
    perf_covar <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), RMSE=numeric(), Rsquare=numeric(), r=numeric(), MAE=numeric())

    # use the same test/train splits as previously defined w/o covariates but rename to avoid writing over
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

    ## train with ridge regression
    ridge_covar <- train(
      as.formula(paste("behav_resids", "~ .")), data = train.data.covar, method = "glmnet",
      #trControl = ctrl_frozen,
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambda_set_covar)
    )

    ## predict in test data
    predictions_ridge_covar <- ridge_covar %>% predict(test.data.covar)
    perf_covar <- data.frame( method = "ridge", brainvar = 'GFC', behavvar=behavvar, nROIs = length(ROIs), iteration = iter, N=nrow(data),
                              RMSE = RMSE(predictions_ridge_covar, test.data.covar$behav_adj),
                              Rsquare = R2(predictions_ridge_covar, test.data.covar$behav_adj),
                              MAE = MAE(predictions_ridge_covar, test.data.covar$behav_adj),
                              r = cor(predictions_ridge_covar, test.data.covar$behav_adj ) )

    ## haufe transform for coefficients per Tian and Zalesky NI 2021
    predictions_ridge_train_covar <- ridge_covar %>% predict(train.data.covar)
    coefs_haufe_covar <- c()
    N_covar <- nrow(train.data.covar)
    for (r in ROIs){ # loop through all edges
      r_std <- scale(train.data.covar[,paste(r)])
      coefs_haufe_covar <- c( coefs_haufe_covar, sum(r_std * predictions_ridge_train_covar) / N_covar )
    }

    ## save out everything
    outname <- paste0(gsub(" ", "_", gsub(":","_",date())), "_", round(runif(1,100,999),0))
    df_covar <- data.frame(snum=test_snums, prediction_ridge=predictions_ridge_covar)
    save(df_covar, coefs_haufe_covar, ROIs_full, file=paste0(workdir4,"/predictions_",outname,".Rdata"))
    save(perf_covar, file=paste0(workdir3,"/performance_",outname,".Rdata"))
  }
}
