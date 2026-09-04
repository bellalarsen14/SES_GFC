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

# set up loop configurations
  # behavvar_idx / covariate_idx index into behavvar_list (1-6)
  # lambda_row / lambda_covar_row index into lambda_selection (row, then column 2)
  # filter_var1 / filter_var2 are the two variables used to restrict to
  #   complete cases before splitting (identical for both members of a pair)
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
motion <- read.csv(paste0(root,'DBIS_GFC_N769_motion.csv'))
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

  workdir <- paste0(root, "Larsen/Dunedin/Bella_Prediction_Outputs/Null_loop/", cfg$name, "/") # this code automatically creates the working folder
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE) # create working directory folder
  
  lambda_set <- lambda_selection[cfg$lambda_row, 2] # from pre-run tuning
  lambda_set_covar <- lambda_selection[cfg$lambda_covar_row, 2]

  # pair-specific complete-case filter, this time within SES level
  behav_merged <- behav_merged_full %>%
    filter(!is.na(.data[[cfg$filter_var1]]) & !is.na(.data[[cfg$filter_var2]]))

  data <- dplyr::left_join(braindat_scaled[,c("snum",ROIs_full)], behav_merged[, c("snum","sex","AverageFD",behavvar)], by="snum")
  data <- data[complete.cases(data),]

  data_covar <- dplyr::left_join(braindat_scaled[,c("snum",ROIs_full)], behav_merged[, c("snum","sex","AverageFD",behavvar,covariate_1)], by="snum")
  data_covar <- data_covar[complete.cases(data_covar),]
  
  # Permutation: shuffle neighborhood deprivation scores before testing prediction
  set.seed(12345) 
  n_perm <- 1000

  r_null       <- numeric(n_perm)
  r_null_covar <- numeric(n_perm)
  
  # -------------------------------------------------------------------------
  # Run permutation
  # -------------------------------------------------------------------------
  
  for (iter in 1:n_perm){
  
    print(paste("iteration", iter, "of", n_perm, "-", cfg$name))
    # first, shuffle SES score in the entire dataset. this will be repeated with different shuffling every iteration.
    data_perm <- data
    data_perm[[behavvar]] <- sample(data[[behavvar]]) # this will re-shuffle with each iteration
    # then create test and train data, this time allowing it to differ between each run, since we are creating 1,000 independent draws
    training.samples <- data[, paste(behavvar)] %>% createDataPartition(p = 0.9, list = FALSE)
    train.data.perm <- data_perm[training.samples, ]
    test.data.perm  <- data_perm[-training.samples, ]
    
    train_snums <- train.data.perm$snum
    test_snums <- test.data.perm$snum
    
    # use the same test/train splits as previously defined w/o covariates but rename to avoid writing over
    data_covar_perm <- data_covar
    data_covar_perm[[behavvar]] <- data_perm[[behavvar]] # keep both models' shuffled outcome identical, since data/data_covar share the same rows
    
    train.data.covar.perm <- data_covar_perm[training.samples, ] # use the same paired 90/10 train/test splits for comparability
    test.data.covar.perm  <- data_covar_perm[-training.samples, ]

    # regress sex and motion from training set
    lm <- lm(train.data.perm[,paste(behavvar)] ~ train.data.perm$sex + train.data.perm$AverageFD)
    train.data.perm$behav_resids <- scale(lm$residuals) 
    # adjust using same parameters in test set
    coefs <- lm$coefficients # 1=intercept, 2=sex, 3=averageFD
    test_fitted_covars <- coefs[1] + coefs[2]*(as.numeric(test.data.perm$sex)-1) +
      coefs[3]*test.data.perm$AverageFD
    test.data.perm$behav_adj <- scale( test.data.perm[,paste(behavvar)] - test_fitted_covars )

    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    train.data.perm <- train.data.perm[, c("behav_resids", ROIs)]
    test.data.perm <- test.data.perm[, c("behav_adj", ROIs)]

    ## train with ridge regression
    ridge <- train(
      as.formula(paste("behav_resids", "~ .")), data = train.data.perm, method = "glmnet",
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambda_set)
    )
    
    # predict in test data
    predictions_ridge_perm <- ridge %>% predict(test.data.perm)
    # save only a vector of the null prediction accuracy values by correlating the null predictions with the true value
    r_null[iter] <- cor(predictions_ridge_perm, test.data.perm$behav_adj)

    #---------------------------------------------------------------------------
    # add covariate and run again
  
    # data is loaded and shuffled above
    
    # regress sex and motion, and covariate_1, from training set
    lm_covar <- lm(train.data.covar.perm[,paste(behavvar)] ~ train.data.covar.perm$sex + train.data.covar.perm$AverageFD + train.data.covar.perm[,paste(covariate_1)])
    train.data.covar.perm$behav_resids <- scale(lm_covar$residuals) 
    # adjust using same parameters in test set
    coefs_covar <- lm_covar$coefficients # 1=intercept, 2=sex, 3=averageFD, 4=covariate
    test_fitted_covars <- coefs_covar[1] + coefs_covar[2]*(as.numeric(test.data.covar.perm$sex)-1) +
      coefs_covar[3]*test.data.covar.perm$AverageFD + coefs_covar[4]*test.data.covar.perm[,paste(covariate_1)]
    test.data.covar.perm$behav_adj <- scale( test.data.covar.perm[,paste(behavvar)] - test_fitted_covars )

    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    train.data.covar.perm <- train.data.covar.perm[, c("behav_resids", ROIs)]
    test.data.covar.perm <- test.data.covar.perm[, c("behav_adj", ROIs)]

    ## train with ridge regression
    ridge_covar <- train(
      as.formula(paste("behav_resids", "~ .")), data = train.data.covar.perm, method = "glmnet",
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = 0, lambda = lambda_set_covar)
    )
    
    # predict in test data
    predictions_ridge_covar_perm <- ridge_covar %>% predict(test.data.covar.perm)
    # save only a vector of the null prediction accuracy values by correlating the null predictions with the true value
    r_null_covar[iter] <- cor(predictions_ridge_covar_perm, test.data.covar.perm$behav_adj)
    
  }
  # save once per config
  save(r_null, file = paste0(workdir, "r_null_", cfg$name, ".Rdata"))
  save(r_null_covar, file = paste0(workdir, "r_null_covar_", cfg$name, ".Rdata"))
}
