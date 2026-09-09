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


citation("apaTables")
# ggseg not yet available for my version of R (4.5.1 as of 02/26/26)
# install.packages("remotes")
# remotes::install_github("ggseg/ggseg")
# remotes::install_github("ggseg/ggsegGlasser")

# install.packages(
#   c("ggseg", "ggseg3d", "ggsegExtra"),
#   repos = c(
#     "https://ggsegverse.r-universe.dev",
#     "https://cloud.r-project.org"
#   )
# )

library(ggseg)
library(ggsegGlasser)

################### OPTIONS #############################
root <- "/parentfoldername/"
setwd("filepathforyourfigures")

# define the behavioral (SES) variables for analysis
behavvar_list <- data.frame(behavvar = c("ADI311","SESchildhd",
                                         "PH45_AreaDeptot","SESall45",
                                         "neighdep2645_factor","ses_composite"))

# load behavioral dataframe

behavdata <- read.csv(file = paste0(root,"person_level_df_all_ses_comp.csv")) %>% 
  select(-X)

# load functional connectivity data
load(paste0(root, 'DBIS_GFC_N769_incSubcortex.Rdata'))
GFC <- FC_ALL
GFC$id <- as.numeric(sub("sub-","",GFC$id)); names(GFC)[1] <- "snum"

# load in motion - average framewise displacement
motion <- read.csv(paste0(root,'DBIS_GFC_N769_motion.csv'))
motion$snum <- as.numeric(sub("sub-","",motion$id))

behavdata <- join_all(list(
  behavdata,
  motion[,c("snum","AverageFD")]
), by="snum",type="full")

# load Glasser parcellation reference file
glasser_to_cab <- read.csv('Glasser_to_CAB.csv', header = TRUE)

# load intraclass correlation (reliability) files for each edge
con_icc <- read.csv('Glasser_ICCs_FCall_DBIS_GFC.FINAL_incSubcortex_mean0.513.csv')

#-------------------------------------------------------------------------------
# Load in performance and prediction outputs for each variable
#-------------------------------------------------------------------------------

# 1: Adult neighborhood-level SES, no covariates

# Performance statistics
file_list_1 <- list.files(paste0(root,'adult_neigh_full/perf_nocovar'))
file_dir_1 <- paste0(root,'adult_neigh_full/perf_nocovar/')

# create an empty dataframe
pred_alpha_df_1 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                            RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

# loop through all of the files in the folder (100), fill in pred_alpha_df
x <- 1
for (f in file_list_1){
  load(paste0(file_dir_1, f))
  pred_alpha_df_1 <- rbind(pred_alpha_df_1, perf)
  x <- x+1
  print(x)
}

# 1a: load in all iterations: PREDICTION
file_list_1a <- list.files(paste0(root,'adult_neigh_full/pred_nocovar'))
file_dir_1a <- paste0(root,'adult_neigh_full/pred_nocovar/')

# create an empty vector to store Haufe-transformed coefficients
feature_importance_neigh_mat <- matrix(rep(0, 100*8805), nrow=100, ncol=8805)
predicted_values_neigh_list <- vector("list", 100)

variable_order_1 <- pred_alpha_df_1$behavvar

# loop through all folders, extract Haufe-transformed coefficents for feature importance
# and extract predicted values per subject 
x <- 1
for (f in file_list_1a){
  load(paste0(file_dir_1a, f))
  if (variable_order_1[x] == 'neighdep2645_factor'){
    feature_importance_neigh_mat[x,] <- c(coefs_haufe) # per model run iteration, extract Haufe-transformed coefficients per edge
    predicted_values_neigh_list[[f]] <- df
    x <- x+1
 
  }
  print(x)
}

# create a feature importance dataframe with the Haufe-transformed coefficients, excluding rows equal to zero
feature_importance_neigh_df <- data.frame(feature_importance_neigh_mat[(rowSums(feature_importance_neigh_mat) != 0),])

# create a dataframe of predicted values per participant across all model runs
predicted_values_neigh_df <- do.call(rbind, predicted_values_neigh_list)

#-------------------------------------------------------------------------------
# 2: Adult individual-level SES, no covariates

# Performance statistics
file_list_2 <- list.files(paste0(root,'adult_ses_full/perf_nocovar'))
file_dir_2 <- paste0(root,'adult_ses_full/perf_nocovar/')

# create an empty dataframe
pred_alpha_df_2 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

# loop through all of the files in the folder (100), fill in pred_alpha_df
x <- 1
for (f in file_list_2){
  load(paste0(file_dir_2, f))
  pred_alpha_df_2 <- rbind(pred_alpha_df_2, perf)
  x <- x+1
  print(x)
}

# 2a: load in all iterations: PREDICTIONS
file_list_2a <- list.files(paste0(root,'adult_ses_full/pred_nocovar'))
file_dir_2a <- paste0(root,'adult_ses_full/pred_nocovar/')

# create an empty vector to store Haufe-transformed coefficients
feature_importance_ses_mat <- matrix(rep(0, 100*8805), nrow=100, ncol=8805)
predicted_values_ses_list <- vector("list", 100)

variable_order_2 <- pred_alpha_df_2$behavvar

# loop through all folders, extract Haufe-transformed coefficents for feature importance
# and extract predicted values per subject 
x <- 1
for (f in file_list_2a){
  load(paste0(file_dir_2a, f))
  if (variable_order_2[x] == 'ses_composite'){
    feature_importance_ses_mat[x,] <- c(coefs_haufe) # per model run iteration, extract Haufe-transformed coefficients per edge
    predicted_values_ses_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

# create a feature importance dataframe with the Haufe-transformed coefficients, excluding rows equal to zero
feature_importance_ses_df <- data.frame(feature_importance_ses_mat[(rowSums(feature_importance_ses_mat) != 0),])

# create a dataframe of predicted values per participant across all model runs
predicted_values_ses_df <- do.call(rbind, predicted_values_ses_list)

#-------------------------------------------------------------------------------
# 3: Childhood individual-level SES, no covariates

# Performance statistics
file_list_3 <- list.files(paste0(root,'chldhd_ses_full/perf_nocovar'))
file_dir_3 <- paste0(root,'chldhd_ses_full/perf_nocovar/')

# create an empty dataframe
pred_alpha_df_3 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

# loop through all of the files in the folder (100), fill in pred_alpha_df
x <- 1
for (f in file_list_3){
  load(paste0(file_dir_3, f))
  pred_alpha_df_3 <- rbind(pred_alpha_df_3, perf)
  x <- x+1
  print(x)
}

# 3a: load in all iterations: PREDICTIONS
file_list_3a <- list.files(paste0(root,'chldhd_ses_full/pred_nocovar'))
file_dir_3a <- paste0(root,'chldhd_ses_full/pred_nocovar/')

# create an empty vector to store Haufe-transformed coefficients
feature_importance_chldhdses_mat <- matrix(rep(0, 100*8805), nrow=100, ncol=8805)
predicted_values_chldhdses_list <- vector("list", 100)

variable_order_3 <- pred_alpha_df_3$behavvar

# loop through all folders, extract Haufe-transformed coefficents for feature importance
# and extract predicted values per subject 
x <- 1
for (f in file_list_3a){
  load(paste0(file_dir_3a, f))
  if (variable_order_3[x] == 'SESchildhd'){
    feature_importance_chldhdses_mat[x,] <- c(coefs_haufe) # per model run iteration, extract Haufe-transformed coefficients per edge
    predicted_values_chldhdses_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

# create a feature importance dataframe with the Haufe-transformed coefficients, excluding rows equal to zero
feature_importance_chldhdses_df <- data.frame(feature_importance_chldhdses_mat[(rowSums(feature_importance_chldhdses_mat) != 0),])

# create a dataframe of predicted values per participant across all model runs
predicted_values_chldhdses_df <- do.call(rbind, predicted_values_chldhdses_list)
#-------------------------------------------------------------------------------
# 4: Childhood neighborhood-level SES, no covariates

# Performance statistics
file_list_4 <- list.files(paste0(root,'chldhd_neigh_full/perf_nocovar'))
file_dir_4 <- paste0(root,'chldhd_neigh_full/perf_nocovar/')

# create an empty dataframe
pred_alpha_df_4 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

# loop through all of the files in the folder (100), fill in pred_alpha_df
x <- 1
for (f in file_list_4){
  load(paste0(file_dir_4, f))
  pred_alpha_df_4 <- rbind(pred_alpha_df_4, perf)
  x <- x+1
  print(x)
}

# 4a: load in all iterations: PREDICTIONS
file_list_4a <- list.files(paste0(root,'chldhd_neigh_full/pred_nocovar'))
file_dir_4a <- paste0(root,'chldhd_neigh_full/pred_nocovar/')

# create an empty vector to store Haufe-transformed coefficients
feature_importance_chldhdneigh_mat <- matrix(rep(0, 100*8805), nrow=100, ncol=8805)
predicted_values_chldhdneigh_list <- vector("list", 100)

variable_order_4 <- pred_alpha_df_4$behavvar

# loop through all folders, extract Haufe-transformed coefficents for feature importance
# and extract predicted values per subject 
x <- 1
for (f in file_list_4a){
  load(paste0(file_dir_4a, f))
  if (variable_order_4[x] == 'ADI311'){
    feature_importance_chldhdneigh_mat[x,] <- c(coefs_haufe) # per model run iteration, extract Haufe-transformed coefficients per edge
    predicted_values_chldhdneigh_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

# create a feature importance dataframe with the Haufe-transformed coefficients, excluding rows equal to zero
feature_importance_chldhdneigh_df <- data.frame(feature_importance_chldhdneigh_mat[(rowSums(feature_importance_chldhdneigh_mat) != 0),])

# create a dataframe of predicted values per participant across all model runs
predicted_values_chldhdneigh_df <- do.call(rbind, predicted_values_chldhdneigh_list)
#-------------------------------------------------------------------------------
# 5: Age 45 individual-level SES, no covariates

# Performance statistics
file_list_5 <- list.files(paste0(root,'age45_ses_full/perf_nocovar'))
file_dir_5 <- paste0(root,'age45_ses_full/perf_nocovar/')

# create an empty dataframe
pred_alpha_df_5 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

# loop through all of the files in the folder (100), fill in pred_alpha_df
x <- 1
for (f in file_list_5){
  load(paste0(file_dir_5, f))
  pred_alpha_df_5 <- rbind(pred_alpha_df_5, perf)
  x <- x+1
  print(x)
}

# 5a: load in all iterations: PREDICTIONS
file_list_5a <- list.files(paste0(root,'age45_ses_full/pred_nocovar'))
file_dir_5a <- paste0(root,'age45_ses_full/pred_nocovar/')

# create an empty vector to store Haufe-transformed coefficients
feature_importance_SES45_mat <- matrix(rep(0, 100*8805), nrow=100, ncol=8805)
predicted_values_ses_45_list <- vector("list", 100)

variable_order_5 <- pred_alpha_df_5$behavvar

# loop through all folders, extract Haufe-transformed coefficents for feature importance
# and extract predicted values per subject 
x <- 1
for (f in file_list_5a){
  load(paste0(file_dir_5a, f))
  if (variable_order_5[x] == 'SESall45'){
    feature_importance_SES45_mat[x,] <- c(coefs_haufe) # per model run iteration, extract Haufe-transformed coefficients per edge
    predicted_values_ses_45_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

# create a feature importance dataframe with the Haufe-transformed coefficients, excluding rows equal to zero
feature_importance_SES45_df <- data.frame(feature_importance_SES45_mat[(rowSums(feature_importance_SES45_mat) != 0),])

# create a dataframe of predicted values per participant across all model runs
predicted_values_ses_45_df <- do.call(rbind, predicted_values_ses_45_list) 
#-------------------------------------------------------------------------------
# 6: Age 45 neighborhood-level SES, no covariates

# Performance statistics
file_list_6 <- list.files(paste0(root,'age45_neigh_full/perf_nocovar'))
file_dir_6 <- paste0(root,'age45_neigh_full/perf_nocovar/')

# create an empty dataframe
pred_alpha_df_6 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                              RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())

# loop through all of the files in the folder (100), fill in pred_alpha_df
x <- 1
for (f in file_list_6){
  load(paste0(file_dir_6, f))
  pred_alpha_df_6 <- rbind(pred_alpha_df_6, perf)
  x <- x+1
  print(x)
}

# 6a: load in all iterations: PREDICTIONS
file_list_6a <- list.files(paste0(root,'age45_neigh_full/pred_nocovar'))
file_dir_6a <- paste0(root,'age45_neigh_full/pred_nocovar/')

# create an empty vector to store Haufe-transformed coefficients
feature_importance_neigh45_mat <- matrix(rep(0, 100*8805), nrow=100, ncol=8805)
predicted_values_neigh_45_list <- vector("list", 100)

variable_order_6 <- pred_alpha_df_6$behavvar

# loop through all folders, extract Haufe-transformed coefficents for feature importance
# and extract predicted values per subject 
x <- 1
for (f in file_list_6a){
  load(paste0(file_dir_6a, f))
  if (variable_order_6[x] == 'PH45_AreaDeptot'){
    feature_importance_neigh45_mat[x,] <- c(coefs_haufe) # per model run iteration, extract Haufe-transformed coefficients per edge
    predicted_values_neigh_45_list[[f]] <- df
    x <- x+1
    
  }
  print(x)
}

# create a feature importance dataframe with the Haufe-transformed coefficients, excluding rows equal to zero 
feature_importance_neigh45_df <- data.frame(feature_importance_neigh45_mat[(rowSums(feature_importance_neigh45_mat) != 0),])

# create a dataframe of predicted values per participant across all model runs
predicted_values_neigh_45_df <- do.call(rbind, predicted_values_neigh_45_list) 

#-------------------------------------------------------------------------------
# Create final feature importance score with all variables by taking the mean of the edges across 100 model iterations (by default excluding ROIs with coef=0)

fi_df <- data.frame(chldhd_SES_fi = colMeans(feature_importance_chldhdses_df),
                    adulthd_ses_fi = colMeans(feature_importance_ses_df),
                    age45_SES_fi = colMeans(feature_importance_SES45_df),
                    chldhd_neigh_fi = colMeans(feature_importance_chldhdneigh_df),
                    adulthd_neigh_fi = colMeans(feature_importance_neigh_df),
                    age45_neigh_fi = colMeans(feature_importance_neigh45_df) 
                    )

# save feature importance df
write.csv(fi_df, "fi_df.csv", row.names = F)

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
  stat_cor(aes(label = after_stat(r.label)),method = "pearson", size = 6,r.digits = 2)+
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
  stat_cor(aes(label = after_stat(r.label)),method = "pearson",size = 6,r.digits = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  theme(text = element_text(size = 18)) +
  labs(title = "Age 45 neighborhood SES (concurrent with scan)",
       x = "Actual neigh. SES",
       y = "Neigh. SES estimated from age 45 FC")

#-------------------------------------------------------------------------------
# create Figure S3 by combining all of the scatterplots
figure_S3 <- ggarrange(chldhd_ses_cor, adult_ses_cor,
                            age45_ses_cor,chldhd_neigh_cor, adult_neigh_cor, age45_neigh_cor)
# annotate_figure(
#   full_cor_plots,
#   top = text_grob("Time of Assessment", size = 16),
#   left = text_grob("Socioeconomic Measure", rot = 90, size = 16)
# )
ggsave("figure_S3.pdf", plot = full_cor_plots, dpi = 300, width = 23, height = 13, units = "in") 

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
  theme(legend.position = "right")

#-------------------------------------------------------------------------------
# Create feature importance heatmaps for each variable (Figure S6)
#-------------------------------------------------------------------------------

# 1: childhood individual-level SES

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_chldhd_ses <- rep(NA, nrow(con_icc))
con_fi_chldhd_ses[con_icc$x > .75] <- colSums(feature_importance_chldhdses_df)

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
png("chldhd_ses.png", height = 900, width = 800, res = 150)

superheat(sd_adj_empty_chldhd_ses,
          membership.rows = network_labels,
          membership.cols = network_labels,
          grid.hline.col = "white",
          grid.vline.col = 'white',
          heat.pal = c('lightblue', 'blue', 'black', 'red', 'yellow'),
          heat.na.col= 'black',
          heat.pal.values = c(0,.5,1),
          heat.lim = c(-5,5),
          title = 'Childhood cumulative SES (birth - 15)',
          bottom.label.text.size = 0,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          left.label.col = network_colors_12,
          legend.height = 0.075,
          legend.text.size = 10)
dev.off()

#-------------------------------------------------------------------------------
# 2: childhood neighborhood-level SES

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_chldhd_neigh <- rep(NA, nrow(con_icc))
con_fi_chldhd_neigh[con_icc$x > .75] <- colSums(feature_importance_chldhdneigh_df)

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

# multiply by -1 for consistent directionality
neg_empty_chldhd_neigh <- -empty_chldhd_neigh

min(neg_empty_chldhd_neigh, na.rm = T)
max(neg_empty_chldhd_neigh, na.rm = T)

max(abs(neg_empty_chldhd_neigh), na.rm=T)

# divide feature importance scores by the standard deviation for visualization
sd_adj_neg_chldhd_neigh <- neg_empty_chldhd_neigh / sd(neg_empty_chldhd_neigh, na.rm = TRUE)

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
png("chldhd_neigh.png", height = 900, width = 800, res = 150)
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
          bottom.label.text.size = 0,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          left.label.col = network_colors_12,
          legend.height = 0.075,
          legend.text.size = 10)
dev.off()

#-------------------------------------------------------------------------------
# 3: adulthood neighborhood-level SES

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_adlthd_neigh <- rep(NA, nrow(con_icc))
con_fi_adlthd_neigh[con_icc$x > .75] <- colSums(feature_importance_neigh_df)

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

# multiply by -1 for consistent directionality
neg_empty_adlthd_neigh <- -empty_adlthd_neigh

# divide feature importance scores by the standard deviation for visualization
sd_adj_neg_adlthd_neigh <- neg_empty_adlthd_neigh / sd(neg_empty_adlthd_neigh, na.rm = TRUE)

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
png("adult_neigh.png", height = 900, width = 800, res = 150)
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
          bottom.label.text.size = 0,
          bottom.label.col = network_colors_12,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          left.label.col = network_colors_12,
          legend.height = 0.075,
          legend.text.size = 10)
dev.off()

#-------------------------------------------------------------------------------
# 4: adulthood individual-level ses

# load in feature importance scores only for edges that are reliable (ICC > 0.75)
con_fi_adlthd_ses <- rep(NA, nrow(con_icc))
con_fi_adlthd_ses[con_icc$x > .75] <- colSums(feature_importance_ses_df)

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
png("adult_ses.png", height = 900, width = 800, res = 150)
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
          bottom.label.text.size = 0,
          bottom.label.text.angle = 90,
          left.label.text.size = 0,
          legend.height = 0.1,
          legend.text.size = 10,
          left.label.col = network_colors_12,
          bottom.label.col = network_colors_12)
dev.off()
