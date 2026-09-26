library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(apaTables)

base_path <- 'Updated_Runs_CV/'

variable_folders <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
# variable_folders should now hold the 6 variable-level folder names (e.g. "chldhd_neigh_full", etc.)

subfolders <- c("haufe_coef_nocovar", "haufe_coef_covar")

coef_list <- list()

for (folder in variable_folders) {
  for (sub in subfolders) {
    dir_path <- file.path(base_path, folder, sub)
    
    rdata_file <- list.files(dir_path, pattern = "\\.[Rr][Dd]ata$",
                             full.names = TRUE)
    
    if (length(rdata_file) == 0) {
      warning(paste("No .Rdata file found in", dir_path))
      next
    }
    if (length(rdata_file) > 1) {
      warning(paste("Multiple .Rdata files in", dir_path, "- using the first"))
      rdata_file <- rdata_file[1]
    }
    
    env <- new.env()
    load(rdata_file, envir = env)
    
    obj_names <- ls(env)
    if (length(obj_names) != 1) {
      warning(paste("Expected exactly 1 object in", rdata_file,
                    "- found", length(obj_names), ":", paste(obj_names, collapse = ", ")))
    }
    
    coef_vec <- as.numeric(get(obj_names[1], envir = env))
    
    suffix <- ifelse(sub == "haufe_coef_nocovar", "nocovar", "covar")
    col_name <- paste(folder, suffix, sep = "_")
    
    coef_list[[col_name]] <- coef_vec
  }
}

# check number of edges: all should read 8805
print(sapply(coef_list, length))

coef_df <- as.data.frame(coef_list)

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

#-------------------------------------------------------------------------------
# Table S7: correlations of feature importance across levels of analysis
cor_mat <- cor(fi_df_nocovar_CV_update, use = "pairwise.complete.obs")

# multiply neighborhood feature importance by -1 for consistent directionality with individual SES
fi_df_nocovar_CV_update <- fi_df_nocovar_CV_update |> 
  mutate(chldhd_neigh_fi_flip = -1*chldhd_neigh_fi,
         adulthd_neigh_fi_flip = -1*adulthd_neigh_fi)

selected_fi <- fi_df_nocovar_CV_update %>% 
  select(chldhd_SES_fi,chldhd_neigh_fi_flip,adulthd_ses_fi,adulthd_neigh_fi_flip)

apa.cor.table(selected_fi, filename = "Table1_APA.doc", table.number = 1)
