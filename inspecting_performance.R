
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)

parent_dir <- 'enteryourparentfoldername'

variable_folders <- list.dirs(parent_dir, full.names = FALSE, recursive = FALSE)

# create an empty dataframe to store results
pred_alpha_df <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                            RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())
x <- 1

for (var_folder in variable_folders) {
  
  perf_folders <- c(
    nocovar = paste0(parent_dir, var_folder, "/perf_nocovar/"),
    covar   = paste0(parent_dir, var_folder, "/perf_covar_ses/")
  )
  
  for (condition in names(perf_folders)) {
    file_dir <- perf_folders[[condition]]
    file_list <- list.files(file_dir)
    
    for (f in file_list) {
      loaded_names <- load(paste0(file_dir, f)) # returns "perf" or "perf_covar" depending on folder
      df_cur <- get(loaded_names[1])
      df_cur$covariate <- condition # tags "nocovar" or "covar"
      df_cur$variable_folder <- var_folder # tags which of the 6 variable folders this came from
      pred_alpha_df <- rbind(pred_alpha_df, df_cur)
      x <- x + 1
      print(x)
    }
  }
}

mean_pred_accuracy <- pred_alpha_df |> 
  group_by(behavvar,covariate) |> 
  summarise(mean(r),mean(Rsquare))

save(mean_pred_accuracy, file = "mean_pred_accuracy.rdata")
write.csv(mean_pred_accuracy, "model_performance.csv", row.names = FALSE)

save(pred_alpha_df, file = "performance_combined_update.rdata")
write.csv(pred_alpha_df, "performance_combined_update.csv", row.names = FALSE)


#-------------------------------------------------------------------------------
# repeat for the models predicting adult SES, covarying for childhood
parent_dir2 <- 'enterparentdirectory2'

variable_folders2 <- list.dirs(parent_dir2, full.names = FALSE, recursive = FALSE)

# create an empty dataframe to store results
pred_alpha_df2 <- data.frame(method=character(), brainvar=character(), behavvar=character(), nROIs=numeric(), iteration=numeric(), N=numeric(), 
                            RMSE=numeric(), Rsquare=numeric(), MAE=numeric(), r=numeric())
x <- 1

for (var_folder in variable_folders2) {
  
  perf_folders <- c(
    covar   = paste0(parent_dir2, var_folder, "/perf_covar_child/")
  )
  
  for (condition in names(perf_folders)) {
    file_dir <- perf_folders[[condition]]
    file_list <- list.files(file_dir)
    
    for (f in file_list) {
      loaded_names <- load(paste0(file_dir, f)) # returns "perf" or "perf_covar" depending on folder
      df_cur <- get(loaded_names[1])
      df_cur$covariate <- condition # tags "nocovar" or "covar"
      df_cur$variable_folder <- var_folder # tags which of the 6 variable folders this came from
      pred_alpha_df2 <- rbind(pred_alpha_df2, df_cur)
      x <- x + 1
      print(x)
    }
  }
}

mean_pred_accuracy2 <- pred_alpha_df2 |> 
  group_by(behavvar,covariate) |> 
  summarise(mean(r),mean(Rsquare))

save(mean_pred_accuracy2, file = "mean_pred_accuracy_covar_child.rdata")
write.csv(mean_pred_accuracy2, "model_performance_covar_child.csv", row.names = FALSE)

