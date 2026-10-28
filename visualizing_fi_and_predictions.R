# load packages
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggExtra)
library(ggplot2)
library(superheat)
library(ggpubr)
library(apaTables)
library(forcats)
library(ggseg)
library(ggseg.formats)
library(stringr)
library(ggseg)
library(ggsegGlasser)

################### OPTIONS #############################
root <- "/Users/uzg9dz/Library/CloudStorage/Box-Box/"
setwd("/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Figures_Updated")

# define the behavioral (SES) variables for analysis
behavvar_list <- data.frame(behavvar = c("ADI311","SESchildhd",
                                         "PH45_AreaDeptot","SESall45",
                                         "neighdep2645_factor","ses_composite"))

# load behavioral dataframe

behavdata <- read.csv(file = "/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Code/person_level_df_all_ses_comp.csv") %>% 
  select(-X)

# load functional connectivity data
load(paste0(root, 'Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_incSubcortex.Rdata'))
GFC <- FC_ALL
GFC$id <- as.numeric(sub("sub-","",GFC$id)); names(GFC)[1] <- "snum"

# load in motion - average framewise displacement
motion <- read.csv("/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_motion.csv")
motion$snum <- as.numeric(sub("sub-","",motion$id))

behavdata <- join_all(list(
  behavdata,
  motion[,c("snum","AverageFD")]
), by="snum",type="full")

# load Glasser parcellation reference file
glasser_to_cab <- read.csv("/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Code\ from\ Ethan/reference\ files/Glasser_to_CAB.csv")

# load intraclass correlation (reliability) files for each edge
con_icc <- read.csv(paste0(root,'Larsen/Dunedin/PH45\ Functional\ Connectivity/DBIS_GFC_N769_incSubcortex_ICCs.csv'))
#-------------------------------------------------------------------------------
# load in feature importance
#-------------------------------------------------------------------------------

base_path <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded'

variable_folders <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
# variable_folders should now hold the 6 variable-level folder names (e.g. "chldhd_neigh_full", etc.)

subfolders <- c("haufe_coef_nocovar", "haufe_coef_covar")

coef_list <- list()
roi_ref   <- NULL   # edge names from the first file; all others must match

for (folder in variable_folders) {
  for (sub in subfolders) {
    dir_path <- file.path(base_path, folder, sub)
    
    rdata_files <- list.files(dir_path, pattern = "\\.[Rr][Dd]ata$", full.names = TRUE)
    if (length(rdata_files) == 0) {
      warning(paste("No .Rdata file found in", dir_path))
      next
    }
    if (length(rdata_files) > 1) {
      warning(paste("Multiple .Rdata files in", dir_path, "- using the most recent"))
    }
    rdata_file <- rdata_files[which.max(file.mtime(rdata_files))]
    
    env <- new.env()
    load(rdata_file, envir = env)
    
    # pick objects by name, not position: matches coefs_haufe and coefs_haufe_covar
    coef_name <- grep("^coefs_haufe", ls(env), value = TRUE)
    if (length(coef_name) != 1 || !exists("ROIs_full", envir = env, inherits = FALSE)) {
      stop("Expected one coefs_haufe* object and ROIs_full in ", rdata_file)
    }
    
    coef_vec <- as.numeric(env[[coef_name]])
    rois     <- env$ROIs_full
    if (length(coef_vec) != length(rois)) {
      stop("Length mismatch in ", rdata_file, ": ", length(coef_vec), " coefs vs ",
           length(rois), " ROIs")
    }
    names(coef_vec) <- rois
    
    if (is.null(roi_ref)) {
      roi_ref <- rois
    } else if (!identical(rois, roi_ref)) {
      stop("ROIs_full in ", rdata_file, " differs from the first file loaded")
    }
    
    suffix   <- ifelse(sub == "haufe_coef_nocovar", "nocovar", "covar")
    col_name <- paste(folder, suffix, sep = "_")
    coef_list[[col_name]] <- coef_vec
  }
}

print(sapply(coef_list, length))   # all should equal length(roi_ref), e.g. 8805

coef_df <- data.frame(coef_list, row.names = roi_ref, check.names = FALSE)

# edge index into the full 78,210-length lower-triangle vector
edge_idx <- as.integer(sub("^edge", "", roi_ref))

coef_df <- coef_df %>% 
  rename(adulthd_neigh_fi_nocovar = adult_neigh_full_nocovar,
         adulthd_neigh_fi_covar = adult_neigh_full_covar,
         adulthd_ses_fi_nocovar = adult_ses_full_nocovar,
         adulthd_ses_fi_covar = adult_ses_full_covar,
         age45_neigh_fi_nocovar = age45_neigh_full_nocovar,
         age45_neigh_fi_covar = age45_neigh_full_covar,
         age45_SES_fi_nocovar = age45_ses_full_nocovar,
         age45_SES_fi_covar = age45_ses_full_covar,
         chldhd_neigh_fi_nocovar = chldhd_neigh_full_nocovar,
         chldhd_neigh_fi_covar = chldhd_neigh_full_covar,
         chldhd_SES_fi_nocovar = chldhd_ses_full_nocovar,
         chldhd_SES_fi_covar = chldhd_ses_full_covar
  )

fi_df_nocovar_CV_update <- coef_df[, grepl("_nocovar$", names(coef_df)) & !grepl("_covar$", names(coef_df))]
names(fi_df_nocovar_CV_update) <- gsub("_nocovar$", "", names(fi_df_nocovar_CV_update))
fi_df <- fi_df_nocovar_CV_update

write.csv(fi_df_nocovar_CV_update, file = "feature_importance_seed_update.csv")
#-------------------------------------------------------------------------------
# Load in performance and prediction outputs for each variable
#-------------------------------------------------------------------------------
# load in one dataframe per variable, with 100 rows and 10 columns. Each row corresponds to one iteration of the ridge regression.
# Columns are the variables saved at each iteration: method, brainvar, nROIs, iteration #, N, RMSE, Rsquare, MAE, and r.
#-------------------------------------------------------------------------------
# 1: adult neighborhood, no covariates
file_list_1 <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_neigh_full/perf_nocovar')
file_dir_1 <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_neigh_full/perf_nocovar/'
pred_alpha_df_1 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

x <- 1
for (f in file_list_1){
  load(paste0(file_dir_1, f))
  pred_alpha_df_1 <- rbind(pred_alpha_df_1, perf)
  x <- x+1
  print(x)
}


# 1a: load in all iterations: PREDICTION
file_list_1a <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_neigh_full/pred_nocovar')
file_dir_1a <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_neigh_full/pred_nocovar/'

predicted_values_neigh_list <- vector("list", 100)

variable_order_1 <- pred_alpha_df_1$behavvar

x <- 1
for (f in file_list_1a){
  load(paste0(file_dir_1a, f))
  if (variable_order_1[x] == 'neighdep2645_factor'){
    predicted_values_neigh_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}
mean(pred_alpha_df_1$r) 

predicted_values_neigh_df <- do.call(rbind, predicted_values_neigh_list)
#-------------------------------------------------------------------------------
# 2: adult ses, no covariates
file_list_2 <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_ses_full/perf_nocovar')
file_dir_2 <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_ses_full/perf_nocovar/'
pred_alpha_df_2 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

x <- 1
for (f in file_list_2){
  load(paste0(file_dir_2, f))
  pred_alpha_df_2 <- rbind(pred_alpha_df_2, perf)
  x <- x+1
  print(x)
}
mean(pred_alpha_df_2$r) 

# 2a: load in all iterations: PREDICTIONS
file_list_2a <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_ses_full/pred_nocovar')
file_dir_2a <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/adult_ses_full/pred_nocovar/'

predicted_values_ses_list <- vector("list", 100)

variable_order_2 <- pred_alpha_df_2$behavvar

x <- 1
for (f in file_list_2a){
  load(paste0(file_dir_2a, f))
  if (variable_order_2[x] == 'ses_composite'){
    predicted_values_ses_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}
 
predicted_values_ses_df <- do.call(rbind, predicted_values_ses_list)
#-------------------------------------------------------------------------------

# 3: Childhood ses, no covariates
file_list_3 <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_ses_full/perf_nocovar')
file_dir_3 <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_ses_full/perf_nocovar/'
pred_alpha_df_3 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

x <- 1
for (f in file_list_3){
  load(paste0(file_dir_3, f))
  pred_alpha_df_3 <- rbind(pred_alpha_df_3, perf)
  x <- x+1
  print(x)
}
mean(pred_alpha_df_3$r) 

# 3a: load in all iterations: PREDICTIONS
file_list_3a <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_ses_full/pred_nocovar')
file_dir_3a <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_ses_full/pred_nocovar/'

predicted_values_chldhdses_list <- vector("list", 100)

variable_order_3 <- pred_alpha_df_3$behavvar

x <- 1
for (f in file_list_3a){
  load(paste0(file_dir_3a, f))
  if (variable_order_3[x] == 'SESchildhd'){
    predicted_values_chldhdses_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

predicted_values_chldhdses_df <- do.call(rbind, predicted_values_chldhdses_list)
#-------------------------------------------------------------------------------
# 4: Childhood neighborhood, no covariates
file_list_4 <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_neigh_full/perf_nocovar')
file_dir_4 <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_neigh_full/perf_nocovar/'
pred_alpha_df_4 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

x <- 1
for (f in file_list_4){
  load(paste0(file_dir_4, f))
  pred_alpha_df_4 <- rbind(pred_alpha_df_4, perf)
  x <- x+1
  print(x)
}
mean(pred_alpha_df_4$r) 

# 4a: load in all iterations: PREDICTIONS
file_list_4a <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_neigh_full/pred_nocovar')
file_dir_4a <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/chldhd_neigh_full/pred_nocovar/'

predicted_values_chldhdneigh_list <- vector("list", 100)

variable_order_4 <- pred_alpha_df_4$behavvar

x <- 1
for (f in file_list_4a){
  load(paste0(file_dir_4a, f))
  if (variable_order_4[x] == 'ADI311'){
    predicted_values_chldhdneigh_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

predicted_values_chldhdneigh_df <- do.call(rbind, predicted_values_chldhdneigh_list)
#-------------------------------------------------------------------------------
# 5: Age 45 SES, no covariates
file_list_5 <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_ses_full/perf_nocovar')
file_dir_5 <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_ses_full/perf_nocovar/'
pred_alpha_df_5 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

x <- 1
for (f in file_list_5){
  load(paste0(file_dir_5, f))
  pred_alpha_df_5 <- rbind(pred_alpha_df_5, perf)
  x <- x+1
  print(x)
}
mean(pred_alpha_df_5$r) 

# 5a: load in all iterations: PREDICTIONS
file_list_5a <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_ses_full/pred_nocovar')
file_dir_5a <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_ses_full/pred_nocovar/'

predicted_values_ses_45_list <- vector("list", 100)

variable_order_5 <- pred_alpha_df_5$behavvar

x <- 1
for (f in file_list_5a){
  load(paste0(file_dir_5a, f))
  if (variable_order_5[x] == 'SESall45'){
    predicted_values_ses_45_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

predicted_values_ses_45_df <- do.call(rbind, predicted_values_ses_45_list) 
#-------------------------------------------------------------------------------
# 6: Age 45 neigh, no covariates
file_list_6 <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_neigh_full/perf_nocovar')
file_dir_6 <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_neigh_full/perf_nocovar/'
pred_alpha_df_6 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

x <- 1
for (f in file_list_6){
  load(paste0(file_dir_6, f))
  pred_alpha_df_6 <- rbind(pred_alpha_df_6, perf)
  x <- x+1
  print(x)
}
mean(pred_alpha_df_6$r) 

# 6a: load in all iterations: PREDICTIONS
file_list_6a <- list.files('/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_neigh_full/pred_nocovar')
file_dir_6a <- '/Users/uzg9dz/Library/CloudStorage/Box-Box/Larsen/Dunedin/Bella_Prediction_Outputs/Updated_Runs_CV_seeded/age45_neigh_full/pred_nocovar/'

predicted_values_neigh_45_list <- vector("list", 100)

variable_order_6 <- pred_alpha_df_6$behavvar

x <- 1
for (f in file_list_6a){
  load(paste0(file_dir_6a, f))
  if (variable_order_6[x] == 'PH45_AreaDeptot'){
    predicted_values_neigh_45_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

predicted_values_neigh_45_df <- do.call(rbind, predicted_values_neigh_45_list) 

#-------------------------------------------------------------------------------
# Create a dataframe of the mean predicted SES value per subject, averaged across model iterations

prediction_df_names <- c("predicted_values_neigh_df","predicted_values_ses_df",
                        "predicted_values_chldhdses_df","predicted_values_chldhdneigh_df",
                        "predicted_values_ses_45_df","predicted_values_neigh_45_df")

predictions_list <- vector("list", length(prediction_df_names))
names(predictions_list) <- prediction_df_names

for (i in seq_along(prediction_df_names)) {
  
  df <- get(prediction_df_names[i])  # pull dataframe by name
  
  predictions_list[[i]] <- df %>%
    group_by(snum) %>%
    summarise(mean_pred = mean(prediction_ridge, na.rm = TRUE),
              .groups = "drop")
}

final_df_predictions <- bind_rows(predictions_list, .id = "source_df")

wide_final_df_predictions <- final_df_predictions %>%
  pivot_wider(
    names_from = source_df,
    values_from = mean_pred
  )

#-------------------------------------------------------------------------------
# To visualize scatterplots of observed and predicted SES scores, merge in real values
behavdata_merge <- behavdata %>% 
  select(snum,sex, neighdep2645_factor,ses_composite,SESchildhd,ADI311,SESall45,PH45_AreaDeptot, AverageFD)

full_pred_actual_df <- inner_join(wide_final_df_predictions,behavdata_merge, by = "snum")

# create scatterplots for each variable
  # note: these correlation values are different than prediction performance r values.
  # These correlations reflect the relationship between observed scores and 
  # averaged predictions, which differ from correlations reported in the main text, which aggregate across participants. 

# 1: adult neighborhood-level SES
adult_neigh_cor <- ggplot(full_pred_actual_df, aes(x=neighdep2645_factor, y=predicted_values_neigh_df)) + 
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)),method = "pearson", size = 6,r.digits = 1)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Adult neighborhood SES (26-45)",
       x = "Actual neigh. SES",
       y = "Neigh. SES estimated from age 45 FC")

# 2: adult individual-level SES
adult_ses_cor <- ggplot(full_pred_actual_df, aes(x=ses_composite, y=predicted_values_ses_df)) + 
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)),method = "pearson",size = 6,r.digits = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Adult individual SES (26-45)",
       x = "Actual SES",
       y = "SES estimated from age 45 FC")

# 3: childhood neighborhood-level SES
chldhd_neigh_cor <- ggplot(full_pred_actual_df, aes(x=ADI311, y=predicted_values_chldhdneigh_df)) + 
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)),method = "pearson",size = 6,r.digits = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Childhood neighborhood SES (3-11)",
       x = "Actual neigh. SES",
       y = "Neigh. SES estimated from age 45 FC")

# 4: childhood individual-level SES
chldhd_ses_cor <- ggplot(full_pred_actual_df, aes(x=SESchildhd, y=predicted_values_chldhdses_df)) + 
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)),method = "pearson",size = 6,r.digits = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Childhood individual SES (birth-15)",
       x = "Actual SES",
       y = "SES estimated from age 45 FC")

# 5: age 45 individual-level SES
age45_ses_cor <- ggplot(full_pred_actual_df, aes(x=SESall45, y=predicted_values_ses_45_df)) + 
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)),method = "pearson",size = 6,r.digits = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Age 45 individual SES (concurrent with scan)",
       x = "Actual SES",
       y = "SES estimated from age 45 FC")

# 5: age 45 neighborhood-level SES
age45_neigh_cor <- ggplot(full_pred_actual_df, aes(x=PH45_AreaDeptot, y=predicted_values_neigh_45_df)) + 
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)),method = "pearson",size = 6,r.digits = 1)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Age 45 neighborhood SES (concurrent with scan)",
       x = "Actual neigh. SES",
       y = "Neigh. SES estimated from age 45 FC")

#-------------------------------------------------------------------------------
# create Figure S3 by combining all of the scatterplots
figure_S3 <- ggarrange(chldhd_ses_cor, adult_ses_cor,
                            age45_ses_cor,chldhd_neigh_cor, adult_neigh_cor, age45_neigh_cor)

ggsave("figure_S3.pdf", plot = figure_S3, dpi = 300, width = 23, height = 13, units = "in") 

#-------------------------------------------------------------------------------
# Table S7: correlations of feature importance across levels of analysis
cor_mat <- cor(fi_df, use = "pairwise.complete.obs")


# multiply neighborhood feature importance by -1 for consistent directionality with individual SES
fi_df <- fi_df |> 
  mutate(chldhd_neigh_fi_flip = -1*chldhd_neigh_fi,
         adulthd_neigh_fi_flip = -1*adulthd_neigh_fi)

selected_fi <- fi_df %>% 
  select(chldhd_SES_fi,chldhd_neigh_fi_flip,adulthd_ses_fi,adulthd_neigh_fi_flip)

apa.cor.table(selected_fi, filename = "Table1_APA.doc", table.number = 1)

cor_mat_spearman <- cor(selected_fi, method = "spearman", use = "pairwise.complete.obs")
round(cor_mat, 2)

#-------------------------------------------------------------------------------

### Visualizations: projecting results onto the connectome

# adding network labels
cab_numbered <- glasser_to_cab$CAB_Network
cab_numbered[cab_numbered == "Visual"] <- 1
cab_numbered[cab_numbered == "Visual2"] <- 2
cab_numbered[cab_numbered == "Somatomotor"] <- 3
cab_numbered[cab_numbered == "Cingulo-Opercular"] <- 4
cab_numbered[cab_numbered == "Language"] <- 5
cab_numbered[cab_numbered == "Frontoparietal"] <- 6
cab_numbered[cab_numbered == "Auditory"] <- 7
cab_numbered[cab_numbered == "Posterior-Multimodal"] <- 8
cab_numbered[cab_numbered == "Default"] <- 9
cab_numbered[cab_numbered == "Dorsal-Attention"] <- 10
cab_numbered[cab_numbered == "Orbito-Affective"] <- 11
cab_numbered[cab_numbered == "Ventral-Multimodal"] <- 12
cab_numbered <- as.numeric(cab_numbered)
cab_numbered <- cab_numbered

network_labels <- glasser_to_cab$CAB_Network
network_labels[is.na(network_labels)] <- "Subcortex"

parcellations <- glasser_to_cab$Glasser_Parcel


# labels with subcortex
network_labels_subcort <-  glasser_to_cab$CAB_Network
network_labels_subcort[is.na(network_labels_subcort)] <- glasser_to_cab$Glasser_Parcel[is.na(network_labels_subcort)]
network_labels_subcort[grepl("Cerebellum", glasser_to_cab$Glasser_Parcel)] <- glasser_to_cab$Glasser_Parcel[grepl("Cerebellum", glasser_to_cab$Glasser_Parcel)]


# Define a color scheme
network_colors_named <- c(
  "Language"              = "#E377C2",
  "Orbito-Affective"      = "#FFBB78",
  "Posterior-Multimodal"  = "#D62728",
  "Somatomotor"           = "#17BECF",
  "Ventral-Multimodal"    = "#FFF300",
  "Visual1"               = "#1F77B4",
  "Dorsal-Attention"      = "#BCBD22", 
  "Frontoparietal"        = "#9467BD",
  "Subcortex"             = "#7F7F7F",
  "Auditory"              = "#AEC7E8",
  "Visual2"               = "#2CA02C",
  "Cingulo-Opercular"     = "#8C564B",
  "Default"               = "#FF7F0E"
)

network_colors_12 <- c("#E377C2", "#FFBB78", "#D62728","#17BECF", "#7F7F7F", "#FFF300",
                       "#1F77B4","#BCBD22","#9467BD","#AEC7E8",
                       "#8C564B","#FF7F0E", "#2CA02C"
)

df_colors <- data.frame(
  network = names(network_colors_named),
  x = 1,
  y = seq_along(network_colors_named)
)

# generate a legend to use with figures
ggplot(df_colors, aes(x = x, y = y, fill = network)) +
  geom_tile() +
  scale_fill_manual(values = network_colors_named) +
  theme_void() +
  theme(legend.position = "top")
#-------------------------------------------------------------------------------
# percentage of positive feature importance scores
fi_df_long <- fi_df |>
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "feature_importance"
  )

fi_df_long <- fi_df_long |> 
  mutate(sign = case_when(feature_importance > 0 ~ "positive",
                          feature_importance < 0 ~ "negative"))
fi_df_long |>
  group_by(variable) |> 
  count(sign) |> 
  mutate(perc = n/sum(n)) |> 
  filter(sign == "positive") |> 
  filter(variable == "adulthd_neigh_fi_flip" |
         variable == "adulthd_ses_fi" |
         variable == "chldhd_SES_fi" |
         variable == "chldhd_neigh_fi_flip")

fi_df_long |>
  count(sign) 

#-------------------------------------------------------------------------------
# Create feature importance heatmaps for each variable (Figure S6)
#-------------------------------------------------------------------------------
# set arguments for the heat pallate
pal  <- c("blue", "#9999FF", "grey96", "#FF9999", "red")
vals <- c(0, 0.35, 0.5, 0.65, 1)

# 1: childhood individual-level SES

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_chldhd_ses <- rep(NA, nrow(con_icc))
con_fi_chldhd_ses[edge_idx] <- fi_df$chldhd_SES_fi

# create an empty matrix (ROI-by-ROI) to save feature importance values
empty_chldhd_ses <- matrix(rep(NA*396*396), nrow=396, ncol=396)

# name rows and columns with Networks
colnames(empty_chldhd_ses) <- network_labels
rownames(empty_chldhd_ses) <- network_labels

# create a parcel-level matrix, where names are Glasser parcellations
glasser_empty_chldhd_ses <- empty_chldhd_ses
colnames(glasser_empty_chldhd_ses) <- parcellations
rownames(glasser_empty_chldhd_ses) <- parcellations

# add fi values to matrices
glasser_empty_chldhd_ses[lower.tri(glasser_empty_chldhd_ses)] <- con_fi_chldhd_ses

empty_chldhd_ses[lower.tri(empty_chldhd_ses)] <- con_fi_chldhd_ses

# flip lower to upper triangles to create symmetrical correlation matrix
t_empty_chldhd_ses <- t(empty_chldhd_ses)[-nrow(empty_chldhd_ses),]
empty_chldhd_ses[upper.tri(empty_chldhd_ses)] <- t_empty_chldhd_ses[upper.tri(t_empty_chldhd_ses, diag = FALSE)]
diag(empty_chldhd_ses) <- NA

min(empty_chldhd_ses, na.rm = T)
max(empty_chldhd_ses, na.rm = T)

# divide feature importance scores by the standard deviation for visualization
sd_adj_empty_chldhd_ses <- empty_chldhd_ses / sd(empty_chldhd_ses, na.rm = TRUE)

min(sd_adj_empty_chldhd_ses, na.rm = T)
max(sd_adj_empty_chldhd_ses, na.rm = T)

superheat(sd_adj_empty_chldhd_ses,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          heat.pal = c('lightblue', 'blue', 'black', 'red', 'yellow'),
          heat.na.col= 'black',
          heat.pal.values = c(0,.5,1),
          heat.lim = c(-5,5),
          title = 'Childhood Individual SES (birth - 15)',
          bottom.label.text.size = 3,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          left.label.col = network_colors_12,
          legend.height = 0.075,
          legend.text.size = 10)

# save a version without labels
tiff("child_ses_heat_sd.tiff", width = 7, height = 7, units = "in", res = 600, compression = "lzw")

superheat(sd_adj_empty_chldhd_ses,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          grid.hline.size = 1.5,
          grid.vline.size = 1.5,
          heat.pal = pal,
          heat.na.col= 'grey97',
          heat.pal.values = vals,
          heat.lim = c(-4,4),
          title = 'Childhood cumulative SES (birth - 15)',
          bottom.label.text.size = 0,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          left.label.col = network_colors_12,
          legend.height = 0,
          legend.text.size = 0,
          bottom.label.size = 0.04,
          left.label.size   = 0.04)
dev.off()

# legend 815, 130
cb <- data.frame(x = seq(-4, 4, length.out = 300), y = 1)

bar <- ggplot(cb, aes(x, y, fill = x)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal, values = vals, limits = c(-4,4), guide = "none") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-4", "-2", "0", "2", "4"),
                     expand = c(0, 0)) +
  labs(x = "SD-scaled Feature importance", y = NULL) +
  theme_minimal(base_size = 20) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(5, 80, 5, 80))   # narrows the bar

#-------------------------------------------------------------------------------
# 2: childhood neighborhood-level SES

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_chldhd_neigh <- rep(NA, nrow(con_icc))
con_fi_chldhd_neigh[edge_idx] <- fi_df$chldhd_neigh_fi_flip

# create an empty matrix (ROI-by-ROI) to save feature importance values
empty_chldhd_neigh <- matrix(rep(NA*396*396), nrow=396, ncol=396)

# name rows and columns with Networks
colnames(empty_chldhd_neigh) <- network_labels
rownames(empty_chldhd_neigh) <- network_labels

# create a parcel-level matrix, where names are Glasser parcellations
glasser_empty_chldhd_neigh <- empty_chldhd_neigh
colnames(glasser_empty_chldhd_neigh) <- parcellations
rownames(glasser_empty_chldhd_neigh) <- parcellations

# add fi values to matrices
glasser_empty_chldhd_neigh[lower.tri(glasser_empty_chldhd_neigh)] <- con_fi_chldhd_neigh

empty_chldhd_neigh[lower.tri(empty_chldhd_neigh)] <- con_fi_chldhd_neigh

# flip lower to upper triangles to create symmetrical correlation matrix
t_empty_chldhd_neigh <- t(empty_chldhd_neigh)[-nrow(empty_chldhd_neigh),]
empty_chldhd_neigh[upper.tri(empty_chldhd_neigh)] <- t_empty_chldhd_neigh[upper.tri(t_empty_chldhd_neigh, diag = FALSE)]
diag(empty_chldhd_neigh) <- NA

min(empty_chldhd_neigh, na.rm = T)
max(empty_chldhd_neigh, na.rm = T)

# divide feature importance scores by the standard deviation for visualization
sd_adj_neg_chldhd_neigh <- empty_chldhd_neigh / sd(empty_chldhd_neigh, na.rm = TRUE)

min(sd_adj_neg_chldhd_neigh, na.rm = T)
max(sd_adj_neg_chldhd_neigh, na.rm = T)

superheat(sd_adj_neg_chldhd_neigh,
          scale = F,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          heat.pal = c('lightblue', 'blue', 'black', 'red', 'yellow'),
          heat.na.col= 'black',
          heat.pal.values = c(0,.5,1),
          heat.lim = c(-5,5),
          title = 'Childhood cumulative neighborhood SES (3 - 11)',
          bottom.label.text.size = 3,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          left.label.col = network_colors_12,
          legend.height = 0.075,
          legend.text.size = 10)

# save a version without labels
tiff("child_neigh_heat_sd.tiff", width = 7, height = 7, units = "in", res = 600, compression = "lzw")
superheat(sd_adj_neg_chldhd_neigh,
          scale = F,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          grid.hline.size = 1.5,
          grid.vline.size = 1.5,
          heat.pal = pal,
          heat.na.col= 'grey97',
          heat.pal.values = vals,
          heat.lim = c(-4.47,4.47),
          title = 'Childhood cumulative neighborhood SES (3 - 11)',
          bottom.label.text.size = 0,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          left.label.col = network_colors_12,
          legend.height = 0,
          legend.text.size = 0,
          bottom.label.size = 0.04,
          left.label.size   = 0.04)
  dev.off()
  
  # legend 815, 130
  cb <- data.frame(x = seq(-4.47, 4.47, length.out = 300), y = 1)
  
  bar <- ggplot(cb, aes(x, y, fill = x)) +
    geom_raster() +
    scale_fill_gradientn(colours = pal, values = vals, limits = c(-4.47,4.47), guide = "none") +
    scale_x_continuous(breaks = c(-4.47, -2.24, 0, 2.24, 4.47),
                       labels = c("-4.47", "-2.24", "0", "2.24", "4.47"),
                       expand = c(0, 0)) +
    labs(x = "SD-scaled Feature importance", y = NULL) +
    theme_minimal(base_size = 20) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          plot.margin = margin(5, 80, 5, 80))   # narrows the bar
  

#-------------------------------------------------------------------------------
# 3: adulthood neighborhood-level SES

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_adlthd_neigh <- rep(NA, nrow(con_icc))
con_fi_adlthd_neigh[edge_idx] <- fi_df$adulthd_neigh_fi_flip

# create an empty matrix (ROI-by-ROI) to save feature importance values
empty_adlthd_neigh <- matrix(rep(NA*396*396), nrow=396, ncol=396)

# name rows and columns with Networks
colnames(empty_adlthd_neigh) <- network_labels
rownames(empty_adlthd_neigh) <- network_labels

# make parcel level matrix
glasser_empty_adlthd_neigh <- empty_adlthd_neigh
colnames(glasser_empty_adlthd_neigh) <- parcellations
rownames(glasser_empty_adlthd_neigh) <- parcellations

# add fi values to matrices
glasser_empty_adlthd_neigh[lower.tri(glasser_empty_adlthd_neigh)] <- con_fi_adlthd_neigh

empty_adlthd_neigh[lower.tri(empty_adlthd_neigh)] <- con_fi_adlthd_neigh

# flip lower to upper triangles
t_empty_adlthd_neigh <- t(empty_adlthd_neigh)[-nrow(empty_adlthd_neigh),]
empty_adlthd_neigh[upper.tri(empty_adlthd_neigh)] <- t_empty_adlthd_neigh[upper.tri(t_empty_adlthd_neigh, diag = FALSE)]
diag(empty_adlthd_neigh) <- NA

min(empty_adlthd_neigh, na.rm = T)
max(empty_adlthd_neigh, na.rm = T)

# divide feature importance scores by the standard deviation for visualization
sd_adj_neg_adlthd_neigh <- empty_adlthd_neigh / sd(empty_adlthd_neigh, na.rm = TRUE)

min(sd_adj_neg_adlthd_neigh, na.rm = T)
max(sd_adj_neg_adlthd_neigh, na.rm = T)

superheat(sd_adj_neg_adlthd_neigh,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          heat.pal = c('lightblue', 'blue', 'black', 'red', 'yellow'),
          heat.na.col= 'black',
          heat.pal.values = c(0,.5,1),
          heat.lim = c(-5,5),
          title = 'Adulthood cumulative neighborhood SES (26 - 46)',
          bottom.label.text.size = 3,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          left.label.col = network_colors_12,
          legend.height = 0.075,
          legend.text.size = 10)

# save a version without labels
tiff("adult_neigh_heat_sd.tiff", width = 7, height = 7, units = "in", res = 600, compression = "lzw")
superheat(sd_adj_neg_adlthd_neigh,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          grid.hline.size = 1.5,
          grid.vline.size = 1.5,
          heat.pal = pal,
          heat.na.col= 'grey97',
          heat.pal.values = vals,
          heat.lim = c(-3.47,3.47),
          title = 'Adulthood cumulative neighborhood SES (26 - 46)',
          bottom.label.text.size = 0,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          left.label.col = network_colors_12,
          legend.height = 0,
          legend.text.size = 0,
          bottom.label.size = 0.04,
          left.label.size   = 0.04)
dev.off()

# legend 815, 130
cb <- data.frame(x = seq(-3.47, 3.47, length.out = 300), y = 1)

bar <- ggplot(cb, aes(x, y, fill = x)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal, values = vals, limits = c(-3.47, 3.47), guide = "none") +
  scale_x_continuous(breaks = c(-3.47, -1.74, 0, 1.74, 3.47),
                     labels = c("-3.47", "-1.74", "0", "1.74", "3.47"),
                     expand = c(0, 0)) +
  labs(x = "SD-scaled Feature importance", y = NULL) +
  theme_minimal(base_size = 20) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(5, 80, 5, 80))   # narrows the bar

#-------------------------------------------------------------------------------
# 4: adulthood individual-level ses

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_adlthd_ses <- rep(NA, nrow(con_icc))
con_fi_adlthd_ses[edge_idx] <- fi_df$adulthd_ses_fi

# create an empty matrix (ROI-by-ROI) to save feature importance values
empty_adlthd_ses <- matrix(rep(NA*396*396), nrow=396, ncol=396)

# name rows and columns with Networks
colnames(empty_adlthd_ses) <- network_labels
rownames(empty_adlthd_ses) <- network_labels

# make parcel level matrix
glasser_adlthd_ses <- empty_adlthd_ses
colnames(glasser_adlthd_ses) <- parcellations
rownames(glasser_adlthd_ses) <- parcellations

# add fi values to matrices
glasser_adlthd_ses[lower.tri(glasser_adlthd_ses)] <- con_fi_adlthd_ses
empty_adlthd_ses[lower.tri(empty_adlthd_ses)] <- con_fi_adlthd_ses

# flip lower to upper triangles
t_empty_adlthd_ses <- t(empty_adlthd_ses)[-nrow(empty_adlthd_ses),]
empty_adlthd_ses[upper.tri(empty_adlthd_ses)] <- t_empty_adlthd_ses[upper.tri(t_empty_adlthd_ses, diag = FALSE)]
diag(empty_adlthd_ses) <- NA

min(empty_adlthd_ses, na.rm = T)
max(empty_adlthd_ses, na.rm = T)

# divide feature importance scores by the standard deviation for visualization
sd_adj_empty_adlthd_ses <- empty_adlthd_ses / sd(empty_adlthd_ses, na.rm = TRUE)

min(sd_adj_empty_adlthd_ses, na.rm = T)
max(sd_adj_empty_adlthd_ses, na.rm = T)

superheat(sd_adj_empty_adlthd_ses,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          heat.pal = c('lightblue', 'blue', 'black', 'red', 'yellow'),
          heat.na.col= 'black',
          heat.pal.values = c(0,.5,1),
          heat.lim = c(-5,5),
          title = 'Adulthood cumulative SES (26 - 46)',
          bottom.label.text.size = 3,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          legend.height = 0.1,
          legend.text.size = 10,
          left.label.col = network_colors_12,
          bottom.label.col = network_colors_12)

# save a version without labels
tiff("adult_ses_heat_sd.tiff", width = 7, height = 7, units = "in", res = 600, compression = "lzw")
superheat(sd_adj_empty_adlthd_ses,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          grid.hline.size = 1.5,
          grid.vline.size = 1.5,
          heat.pal = pal,
          heat.na.col= 'grey97',
          heat.pal.values = vals,
          heat.lim = c(-4.9,4.9),
          title = 'Adulthood cumulative SES (26 - 46)',
          bottom.label.text.size = 0,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          legend.height = 0,
          legend.text.size = 0,
          left.label.col = network_colors_12,
          bottom.label.col = network_colors_12,
          bottom.label.size = 0.04,
          left.label.size   = 0.04)
dev.off()

# legend 815, 130
cb <- data.frame(x = seq(-4.9,4.9, length.out = 300), y = 1)

bar <- ggplot(cb, aes(x, y, fill = x)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal, values = vals, limits = c(-4.9,4.9), guide = "none") +
  scale_x_continuous(breaks = c(-4.9, -2.45, 0, 2.45, 4.9),
                     labels = c("-4.9", "-2.45", "0", "2.45", "4.9"),
                     expand = c(0, 0)) +
  labs(x = "SD-scaled Feature importance", y = NULL) +
  theme_minimal(base_size = 20) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(5, 80, 5, 80))   # narrows the bar
#-------------------------------------------------------------------------------
# legend, for all
cb <- data.frame(x = seq(lim[1], lim[2], length.out = 300), y = 1)
bar <- ggplot(cb, aes(x, y, fill = x)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal, values = vals, limits = lim, guide = "none") +
  scale_x_continuous(breaks = c(-0.14, -0.07, 0, 0.07, 0.14),
                     labels = c("-0.14", "-0.07", "0", "0.07", "0.14"),
                     expand = c(0, 0)) +
  labs(x = "SD-scaled Feature importance", y = NULL) +
  theme_minimal(base_size = 20) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(5, 80, 5, 80))   # narrows the bar

#-------------------------------------------------------------------------------
# nodal level
#-------------------------------------------------------------------------------
# childhood ses

# make parcel level matrix
sd_adj_glasser_chldhd_ses <- glasser_empty_chldhd_ses / sd(glasser_empty_chldhd_ses, na.rm = TRUE)

sd_adj_glasser_chldhd_ses[is.na(sd_adj_glasser_chldhd_ses)] <- 0  # make NA's 0

sd_nodal_fi_chldhd_ses <-  data.frame(Glasser_Parcel = colnames(sd_adj_glasser_chldhd_ses),
                                      chldhd_ses_fi = colMeans(sd_adj_glasser_chldhd_ses, na.rm=T))

sd_nodal_fi_chldhd_ses <- full_join(sd_nodal_fi_chldhd_ses,glasser_to_cab, by = "Glasser_Parcel")

#-------------------------------------------------------------------------------
# childhood neighborhood

# make parcel level matrix
sd_adj_glasser_chldhd_neigh <- glasser_empty_chldhd_neigh / sd(glasser_empty_chldhd_neigh, na.rm = TRUE)

sd_adj_glasser_chldhd_neigh[is.na(sd_adj_glasser_chldhd_neigh)] <- 0  # make NA's 0

sd_nodal_fi_chldhd_neigh <-  data.frame(Glasser_Parcel = colnames(sd_adj_glasser_chldhd_neigh),
                                      chldhd_neigh_fi = colMeans(sd_adj_glasser_chldhd_neigh, na.rm=T))

sd_nodal_fi_chldhd_neigh <- full_join(sd_nodal_fi_chldhd_neigh,glasser_to_cab, by = "Glasser_Parcel")

#-------------------------------------------------------------------------------
# adult ses

# make parcel level matrix
sd_adj_glasser_adult_ses <- glasser_adlthd_ses / sd(glasser_adlthd_ses, na.rm = TRUE)

sd_adj_glasser_adult_ses[is.na(sd_adj_glasser_adult_ses)] <- 0  # make NA's 0

sd_nodal_fi_adult_ses <-  data.frame(Glasser_Parcel = colnames(sd_adj_glasser_adult_ses),
                                        adult_ses_fi = colMeans(sd_adj_glasser_adult_ses, na.rm=T))

sd_nodal_fi_adult_ses <- full_join(sd_nodal_fi_adult_ses,glasser_to_cab, by = "Glasser_Parcel")

#-------------------------------------------------------------------------------
# adult neighborhood

# make parcel level matrix
sd_adj_glasser_adult_neigh <- glasser_empty_adlthd_neigh / sd(glasser_empty_adlthd_neigh, na.rm = TRUE)

sd_adj_glasser_adult_neigh[is.na(sd_adj_glasser_adult_neigh)] <- 0  # make NA's 0

sd_nodal_fi_adult_neigh <-  data.frame(Glasser_Parcel = colnames(sd_adj_glasser_adult_neigh),
                                     adult_neigh_fi = colMeans(sd_adj_glasser_adult_neigh, na.rm=T))

sd_nodal_fi_adult_neigh <- full_join(sd_nodal_fi_adult_neigh,glasser_to_cab, by = "Glasser_Parcel")


