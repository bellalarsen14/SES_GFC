# load packages
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(broom)

# estimate statistical significance by comparing null results to observed results

root <- "/parentfoldername/"
num_perm <- 1000

# load in true model results, the average prediction accuracy value per variable
load("mean_pred_accuracy.rdata")

#-------------------------------------------------------------------------------
# providing an example of one SES variable; same code is repeated to calculate the p-value of each variable

# CHILDHOOD NEIGHBORHOOD-LEVEL SES

# load in NULL iterations from permutation
load(paste0(root, 'chldhd_neigh/r_null_chldhd_neigh.Rdata'))
r_null_chldhd_neigh_nocovar <- r_null

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_chldhd_neigh_nocovar <- pull(mean_pred_accuracy[2,3])
p_chldhd_neigh_nocovar <- length(which(r_null_chldhd_neigh_nocovar>=true_mean_r_chldhd_neigh_nocovar))/num_perm
null_mean_r_chldhd_neigh_nocovar <- mean(r_null_chldhd_neigh_nocovar)
# calculate percentile
diff_mean_chldhd_neigh_nocovar <- true_mean_r_chldhd_neigh_nocovar - null_mean_r_chldhd_neigh_nocovar
percentile_0 <- pnorm(diff_mean_chldhd_neigh_nocovar, mean = null_mean_r_chldhd_neigh_nocovar, sd = sd(r_null_chldhd_neigh_nocovar)) * 100

# save p-value
p_value_0 <- data.frame(
  behavvar = c("ADI311"),
  covariate = c("nocovar"),
  p_value = p_chldhd_neigh_nocovar,
  percentile = percentile_0
)

# COVARYING FOR CHILDHOOD INDIVIDUAL-LEVEL SES
# load in NULL iterations
load(paste0(root, 'chldhd_neigh/r_null_covar_chldhd_neigh.Rdata'))
r_null_chldhd_neigh_covar <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_chldhd_neigh_covar <- pull(mean_pred_accuracy[1,3])
p_chldhd_neigh_covar <- length(which(r_null_chldhd_neigh_covar>=true_mean_r_chldhd_neigh_covar))/num_perm
null_mean_r_chldhd_neigh_covar <- mean(r_null_chldhd_neigh_covar)
# calculate percentile
diff_mean_chldhd_neigh_covar <- true_mean_r_chldhd_neigh_covar - null_mean_r_chldhd_neigh_covar
percentile_1 <- pnorm(diff_mean_chldhd_neigh_covar, mean = null_mean_r_chldhd_neigh_covar, sd = sd(r_null_chldhd_neigh_covar)) * 100

# save p-value
p_value_1 <- data.frame(
  behavvar = c("ADI311"),
  covariate = c("covar"),
  p_value = p_chldhd_neigh_covar,
  percentile = percentile_1
)

#-------------------------------------------------------------------------------
# Models predicting adulthood SES, covarying for childhood

  # load in TRUE accuracy
load("mean_pred_accuracy_covar_child.rdata")

#-------------------------------------------------------------------------------
# ADULT NEIGHBORHOOD-LEVEL SES, COVARYING FOR CHILD NEIGHBORHOOD-LEVEL SES

# load in NULL iterations
load("r_null_covar_adult_neigh.Rdata")
r_null_adult_neigh_covar_child <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_adult_neigh_covar_child <- pull(mean_pred_accuracy2[1,3])
p_adult_neigh_covar_child <- length(which(r_null_adult_neigh_covar_child>=true_mean_r_adult_neigh_covar_child))/num_perm
null_mean_r_adult_neigh_covar_child <- mean(r_null_adult_neigh_covar_child)
# calculate percentile
diff_mean_adult_neigh_covar_child <- true_mean_r_adult_neigh_covar_child - null_mean_r_adult_neigh_covar_child
percentile_12 <- pnorm(diff_mean_adult_neigh_covar_child, mean = null_mean_r_adult_neigh_covar_child, sd = sd(r_null_adult_neigh_covar_child)) * 100

# save p-value
p_value_12 <- data.frame(
  behavvar = c("neighdep2645_factor"),
  covariate = c("covar_child"),
  p_value = p_adult_neigh_covar_child,
  percentile = percentile_12
)

#-------------------------------------------------------------------------------
# ADULT INDIVIDUAL-LEVEL SES, COVARYING FOR CHILD INDIVIDUAL-LEVEL SES

# load in NULL iterations
load("r_null_covar_adult_ses.Rdata")
r_null_adult_ses_covar_child <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_adult_ses_covar_child <- pull(mean_pred_accuracy2[2,3])
p_adult_ses_covar_child <- length(which(r_null_adult_ses_covar_child>=true_mean_r_adult_ses_covar_child))/num_perm
null_mean_r_adult_ses_covar_child <- mean(r_null_adult_ses_covar_child)
# calculate percentile
diff_mean_adult_ses_covar_child <- true_mean_r_adult_ses_covar_child - null_mean_r_adult_ses_covar_child
percentile_13 <- pnorm(diff_mean_adult_ses_covar_child, mean = null_mean_r_adult_ses_covar_child, sd = sd(r_null_adult_ses_covar_child)) * 100

# save p-value
p_value_13 <- data.frame(
  behavvar = c("ses_composite"),
  covariate = c("covar_child"),
  p_value = p_adult_ses_covar_child,
  percentile = percentile_13
)

#-------------------------------------------------------------------------------
# merge with true results
mean_pred_accuracy2$covariate <- c("covar_child","covar_child")
mean_pred_full <- rbind(mean_pred_accuracy,mean_pred_accuracy2)
results_full <- full_join(mean_pred_full,p_value_full)
results_full <- results_full |> 
  group_by(covariate) |> 
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))

