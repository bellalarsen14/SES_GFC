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

# load in true model results
load("mean_pred_accuracy.rdata")

#-------------------------------------------------------------------------------
# CHILDHOOD NEIGHBORHOOD-LEVEL SES

# load in NULL iterations
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

p_value_full <- rbind(p_value_0,p_value_1)
#-------------------------------------------------------------------------------
# CHILDHOOD INDIVIDUAL-LEVEL SES

# load in NULL iterations
load(paste0(root, 'chldhd_ses/r_null_chldhd_ses.Rdata'))
r_null_chldhd_ses_nocovar <- r_null

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_chldhd_ses_nocovar <- pull(mean_pred_accuracy[8,3])
p_chldhd_ses_nocovar <- length(which(r_null_chldhd_ses_nocovar>=true_mean_r_chldhd_ses_nocovar))/num_perm
null_mean_r_chldhd_ses_nocovar <- mean(r_null_chldhd_ses_nocovar)
# calculate percentile
diff_mean_chldhd_ses_nocovar <- true_mean_r_chldhd_ses_nocovar - null_mean_r_chldhd_ses_nocovar
percentile_2 <- pnorm(diff_mean_chldhd_ses_nocovar, mean = null_mean_r_chldhd_ses_nocovar, sd = sd(r_null_chldhd_ses_nocovar)) * 100

# save p-value
p_value_2 <- data.frame(
  behavvar = c("SESchildhd"),
  covariate = c("nocovar"),
  p_value = p_chldhd_ses_nocovar,
  percentile = percentile_2
)

# COVARYING FOR CHILDHOOD NEIGHBORHOOD-LEVEL SES
# load in NULL iterations
load(paste0(root, 'chldhd_ses/r_null_covar_chldhd_ses.Rdata'))
r_null_chldhd_ses_covar <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_chldhd_ses_covar <- pull(mean_pred_accuracy[7,3])
p_chldhd_ses_covar <- length(which(r_null_chldhd_ses_covar>=true_mean_r_chldhd_ses_covar))/num_perm
null_mean_r_chldhd_ses_covar <- mean(r_null_chldhd_ses_covar)
# calculate percentile
diff_mean_chldhd_ses_covar <- true_mean_r_chldhd_ses_covar - null_mean_r_chldhd_ses_covar
percentile_3 <- pnorm(diff_mean_chldhd_ses_covar, mean = null_mean_r_chldhd_ses_covar, sd = sd(r_null_chldhd_ses_covar)) * 100

# save p-value
p_value_3 <- data.frame(
  behavvar = c("SESchildhd"),
  covariate = c("covar"),
  p_value = p_chldhd_ses_covar,
  percentile = percentile_3
)

p_value_full <- rbind(p_value_full,p_value_2,p_value_3)
#-------------------------------------------------------------------------------
# ADULT NEIGHBORHOOD-LEVEL SES

# load in NULL iterations
load(paste0(root, 'adult_neigh/r_null_adult_neigh.Rdata'))
r_null_adult_neigh_nocovar <- r_null

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_adult_neigh_nocovar <- pull(mean_pred_accuracy[10,3])
p_adult_neigh_nocovar <- length(which(r_null_adult_neigh_nocovar>=true_mean_r_adult_neigh_nocovar))/num_perm
null_mean_r_adult_neigh_nocovar <- mean(r_null_adult_neigh_nocovar)
# calculate percentile
diff_mean_adult_neigh_nocovar <- true_mean_r_adult_neigh_nocovar - null_mean_r_adult_neigh_nocovar
percentile_4 <- pnorm(diff_mean_adult_neigh_nocovar, mean = null_mean_r_adult_neigh_nocovar, sd = sd(r_null_adult_neigh_nocovar)) * 100

# save p-value
p_value_4 <- data.frame(
  behavvar = c("neighdep2645_factor"),
  covariate = c("nocovar"),
  p_value = p_adult_neigh_nocovar,
  percentile = percentile_4
)

# COVARYING FOR ADULTHOOD INDIVIDUAL-LEVEL SES

# load in NULL iterations
load(paste0(root, 'adult_neigh/r_null_covar_adult_neigh.Rdata'))
r_null_adult_neigh_covar <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_adult_neigh_covar <- pull(mean_pred_accuracy[9,3])
p_adult_neigh_covar <- length(which(r_null_adult_neigh_covar>=true_mean_r_adult_neigh_covar))/num_perm
null_mean_r_adult_neigh_covar <- mean(r_null_adult_neigh_covar)
# calculate percentile
diff_mean_adult_neigh_covar <- true_mean_r_adult_neigh_covar - null_mean_r_adult_neigh_covar
percentile_5 <- pnorm(diff_mean_adult_neigh_covar, mean = null_mean_r_adult_neigh_covar, sd = sd(r_null_adult_neigh_covar)) * 100

# save p-value
p_value_5 <- data.frame(
  behavvar = c("neighdep2645_factor"),
  covariate = c("covar"),
  p_value = p_adult_neigh_covar,
  percentile = percentile_5
)

p_value_full <- rbind(p_value_full,p_value_4,p_value_5)

#-------------------------------------------------------------------------------
# ADULTHOOD INDIVIDUAL-LEVEL SES

# load in NULL iterations
load(paste0(root, 'adult_ses/r_null_adult_ses.Rdata'))
r_null_adult_ses_nocovar <- r_null

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_adult_ses_nocovar <- pull(mean_pred_accuracy[12,3])
p_adult_ses_nocovar <- length(which(r_null_adult_ses_nocovar>=true_mean_r_adult_ses_nocovar))/num_perm
null_mean_r_adult_ses_nocovar <- mean(r_null_adult_ses_nocovar)
# calculate percentile
diff_mean_adult_ses_nocovar <- true_mean_r_adult_ses_nocovar - null_mean_r_adult_ses_nocovar
percentile_6 <- pnorm(diff_mean_adult_ses_nocovar, mean = null_mean_r_adult_ses_nocovar, sd = sd(r_null_adult_ses_nocovar)) * 100

# save p-value
p_value_6 <- data.frame(
  behavvar = c("ses_composite"),
  covariate = c("nocovar"),
  p_value = p_adult_ses_nocovar,
  percentile = percentile_6
)

# COVARYING FOR ADULTHOOD NEIGHBORHOOD-LEVEL SES
# load in NULL iterations
load(paste0(root, 'adult_ses/r_null_covar_adult_ses.Rdata'))
r_null_adult_ses_covar <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_adult_ses_covar <- pull(mean_pred_accuracy[11,3])
p_adult_ses_covar <- length(which(r_null_adult_ses_covar>=true_mean_r_adult_ses_covar))/num_perm
null_mean_r_adult_ses_covar <- mean(r_null_adult_ses_covar)
# calculate percentile
diff_mean_adult_ses_covar <- true_mean_r_adult_ses_covar - null_mean_r_adult_ses_covar
percentile_7 <- pnorm(diff_mean_adult_ses_covar, mean = null_mean_r_adult_ses_covar, sd = sd(r_null_adult_ses_covar)) * 100

# save p-value
p_value_7 <- data.frame(
  behavvar = c("ses_composite"),
  covariate = c("covar"),
  p_value = p_adult_ses_covar,
  percentile = percentile_7
)

p_value_full <- rbind(p_value_full,p_value_6,p_value_7)

#-------------------------------------------------------------------------------
# AGE 45 NEIGHBORHOOD-LEVEL SES

# load in NULL iterations
load(paste0(root, 'age45_neigh/r_null_age45_neigh.Rdata'))
r_null_age45_neigh_nocovar <- r_null

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_age45_neigh_nocovar <- pull(mean_pred_accuracy[4,3])
p_age45_neigh_nocovar <- length(which(r_null_age45_neigh_nocovar>=true_mean_r_age45_neigh_nocovar))/num_perm
null_mean_r_age45_neigh_nocovar <- mean(r_null_age45_neigh_nocovar)
# calculate percentile
diff_mean_age45_neigh_nocovar <- true_mean_r_age45_neigh_nocovar - null_mean_r_age45_neigh_nocovar
percentile_8 <- pnorm(diff_mean_age45_neigh_nocovar, mean = null_mean_r_age45_neigh_nocovar, sd = sd(r_null_age45_neigh_nocovar)) * 100

# save p-value
p_value_8 <- data.frame(
  behavvar = c("PH45_AreaDeptot"),
  covariate = c("nocovar"),
  p_value = p_age45_neigh_nocovar,
  percentile = percentile_8
)

# COVARYING FOR AGE 45 INDIVIDUAL-LEVEL SES
# load in NULL iterations
load(paste0(root, 'age45_neigh/r_null_covar_age45_neigh.Rdata'))
r_null_age45_neigh_covar <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_age45_neigh_covar <- pull(mean_pred_accuracy[3,3])
p_age45_neigh_covar <- length(which(r_null_age45_neigh_covar>=true_mean_r_age45_neigh_covar))/num_perm
null_mean_r_age45_neigh_covar <- mean(r_null_age45_neigh_covar)
# calculate percentile
diff_mean_age45_neigh_covar <- true_mean_r_age45_neigh_covar - null_mean_r_age45_neigh_covar
percentile_9 <- pnorm(diff_mean_age45_neigh_covar, mean = null_mean_r_age45_neigh_covar, sd = sd(r_null_age45_neigh_covar)) * 100

# save p-value
p_value_9 <- data.frame(
  behavvar = c("PH45_AreaDeptot"),
  covariate = c("covar"),
  p_value = p_age45_neigh_covar,
  percentile = percentile_9
)

p_value_full <- rbind(p_value_full,p_value_8,p_value_9)

#-------------------------------------------------------------------------------
# AGE 45 INDIVIDUAL-LEVEL SES

# load in NULL iterations
load(paste0(root, 'age45_ses/r_null_age45_ses.Rdata'))
r_null_age45_ses_nocovar <- r_null

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_age45_ses_nocovar <- pull(mean_pred_accuracy[6,3])
p_age45_ses_nocovar <- length(which(r_null_age45_ses_nocovar>=true_mean_r_age45_ses_nocovar))/num_perm
null_mean_r_age45_ses_nocovar <- mean(r_null_age45_ses_nocovar)
# calculate percentile
diff_mean_age45_ses_nocovar <- true_mean_r_age45_ses_nocovar - null_mean_r_age45_ses_nocovar
percentile_10 <- pnorm(diff_mean_age45_ses_nocovar, mean = null_mean_r_age45_ses_nocovar, sd = sd(r_null_age45_ses_nocovar)) * 100

# save p-value
p_value_10 <- data.frame(
  behavvar = c("SESall45"),
  covariate = c("nocovar"),
  p_value = p_age45_ses_nocovar,
  percentile = percentile_10
)

# COVARYING FOR AGE 45 NEIGHBORHOOD-LEVEL SES
# load in NULL iterations
load(paste0(root, 'age45_ses/r_null_covar_age45_ses.Rdata'))
r_null_age45_ses_covar <- r_null_covar

# calculate p-value: proportion of null model r's greater than or equal to those corresponding to the original mean r
true_mean_r_age45_ses_covar <- pull(mean_pred_accuracy[5,3])
p_age45_ses_covar <- length(which(r_null_age45_ses_covar>=true_mean_r_age45_ses_covar))/num_perm
null_mean_r_age45_ses_covar <- mean(r_null_age45_ses_covar)
# calculate percentile
diff_mean_age45_ses_covar <- true_mean_r_age45_ses_covar - null_mean_r_age45_ses_covar
percentile_11 <- pnorm(diff_mean_age45_ses_covar, mean = null_mean_r_age45_ses_covar, sd = sd(r_null_age45_ses_covar)) * 100

# save p-value
p_value_11 <- data.frame(
  behavvar = c("SESall45"),
  covariate = c("covar"),
  p_value = p_age45_ses_covar,
  percentile = percentile_11
)

p_value_full <- rbind(p_value_full,p_value_10,p_value_11)

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

p_value_full <- rbind(p_value_full,p_value_12,p_value_13)

#-------------------------------------------------------------------------------
# merge with true results
mean_pred_accuracy2$covariate <- c("covar_child","covar_child")
mean_pred_full <- rbind(mean_pred_accuracy,mean_pred_accuracy2)
results_full <- full_join(mean_pred_full,p_value_full)
results_full <- results_full |> 
  group_by(covariate) |> 
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))

# save
write.csv(results_full, "results_full.csv", row.names = FALSE)
