# load packages
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(broom)
library(ggrain)
library(ggsignif)
library(ggpubr)

################### OPTIONS #############################
root <- "/parentfoldername/"

# define the behavioral (SES) variables for analysis
behavvar_list <- data.frame(behavvar = c("ADI311","SESchildhd",
                                         "PH45_AreaDeptot","SESall45",
                                         "neighdep2645_factor","ses_composite"))

# load performance statistics
performance_combined <- read.csv(paste0(root,'performance_combined_update.csv')) 

# rename variables for interpretability
performance_combined <- performance_combined %>%
  mutate(ses_var = case_when(behavvar == "neighdep2645_factor" ~ "Neighborhood Deprivation",
                             behavvar == "PH45_AreaDeptot" ~ "Neighborhood Deprivation",
                             behavvar == "ADI311" ~ "Neighborhood Deprivation",
                             behavvar == "ses_composite" ~ "SES",
                             behavvar == "SESall45" ~ "SES",
                             behavvar == "SESchildhd" ~ "SES")) |> 
  mutate(time_assessed = case_when(behavvar == "neighdep2645_factor" ~ "Adulthood",
                             behavvar == "PH45_AreaDeptot" ~ "Age 45",
                             behavvar == "ADI311" ~ "Childhood",
                             behavvar == "ses_composite" ~ "Adulthood",
                             behavvar == "SESall45" ~ "Age 45",
                             behavvar == "SESchildhd" ~ "Childhood")) |>
  mutate(covariates_yn = case_when(covariate == "nocovar" ~ "N",
                                   covariate == "covar" ~ "Y")) |>
  mutate(ses_var_rename = case_when(ses_var == "Neighborhood Deprivation" ~ "Neighborhood SES",
                                    ses_var == "SES" ~ "Individual SES"),
         time_assessed_rename = case_when(time_assessed == "Childhood" ~ "Childhood (birth-15)",
                                          time_assessed == "Adulthood" ~ "Adulthood (26-45)",
                                          time_assessed == "Age 45" ~ "Age 45 (concurrent with scan)"),
         time_assessed_split = case_when(time_assessed == "Childhood" & ses_var == "Neighborhood Deprivation" ~ "Childhood (3-11)",
                                         time_assessed == "Childhood" & ses_var == "SES" ~ "Childhood (birth-15)",
                                         time_assessed == "Adulthood" ~ "Adulthood (26-45)",
                                         time_assessed == "Age 45" ~ "Age 45 (concurrent with scan)"))

# between-group statistical significance 

results_list <- list()

counter <- 1

for (t in unique(performance_combined$time_assessed)) {
  for (c in unique(performance_combined$covariates_yn)) {
    
    subset_data <- subset(performance_combined, 
                          time_assessed == t & covariates_yn == c)
    
    test <- t.test(r ~ ses_var_rename, data = subset_data)
    
    results_list[[counter]] <- data.frame(
      time_assessed = t,
      covariates_yn  = c,
      t_value    = as.numeric(test$statistic),
      df         = as.numeric(test$parameter),
      p_value    = test$p.value,
      conf_low   = test$conf.int[1],
      conf_high  = test$conf.int[2]
    )
    
    counter <- counter + 1
  }
}

results_df <- do.call(rbind, results_list)

results_df

results_df <- results_df %>%
  group_by(covariates_yn) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

write.csv(results_df, file = "between_group_results_df.csv", row.names = TRUE)
#-------------------------------------------------------------------------------

# Individual ses between-age statistical significance (create loop through all combos)
# test manually

individ_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% filter(covariates_yn == "N" & time_assessed != "Age 45" & ses_var_rename == "Individual SES"))
neigh_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% filter(covariates_yn == "N" & time_assessed != "Age 45" & ses_var_rename == "Neighborhood SES"))


results_list_t_separated <- data.frame(
  ses_var = c("Individual SES","Neighborhood SES"),
  t_value    = c(as.numeric(individ_ses_t$statistic),as.numeric(neigh_ses_t$statistic)),
  df         = c(as.numeric(individ_ses_t$parameter),as.numeric(neigh_ses_t$parameter)),
  p_value    = c(individ_ses_t$p.value,neigh_ses_t$p.value),
  conf_low   = c(individ_ses_t$conf.int[1],neigh_ses_t$conf.int[1]),
  conf_high  = c(individ_ses_t$conf.int[2],neigh_ses_t$conf.int[2])
  )

results_child_vs_adult <- results_list_t_separated %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) 

write.csv(results_child_vs_adult, file = "results_child_vs_adult.csv", row.names = TRUE)

#-------------------------------------------------------------------------------
# age 45 results

individ_45_adult_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% 
                                   filter(covariates_yn == "N" & time_assessed != "Childhood" & ses_var_rename == "Individual SES"))
individ_45_child_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% 
                                   filter(covariates_yn == "N" & time_assessed != "Adulthood" & ses_var_rename == "Individual SES"))
neigh_45_adult_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% 
                                 filter(covariates_yn == "N" & time_assessed != "Childhood" & ses_var_rename == "Neighborhood SES"))
neigh_45_child_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% 
                                 filter(covariates_yn == "N" & time_assessed != "Adulthood" & ses_var_rename == "Neighborhood SES"))


results_list_t_separated_full <- data.frame(
  ses_var = c("Individual SES (childhood vs. adulthood)",
              "Neighborhood SES (childhood vs. adulthood)",
              "Individual SES (age 45 vs. childhood)",
              "Neighborhood SES (age 45 vs. childhood)",
              "Individual SES (age 45 vs. adulthood)",
              "Neighborhood SES (age 45 vs. adulthood)"),
  t_value    = c(as.numeric(individ_ses_t$statistic),as.numeric(neigh_ses_t$statistic),
                 as.numeric(individ_45_child_ses_t$statistic),as.numeric(neigh_45_child_ses_t$statistic),
                 as.numeric(individ_45_adult_ses_t$statistic),as.numeric(neigh_45_adult_ses_t$statistic)),
  df         = c(as.numeric(individ_ses_t$parameter),as.numeric(neigh_ses_t$parameter),
                 as.numeric(individ_45_child_ses_t$parameter),as.numeric(neigh_45_child_ses_t$parameter),
                 as.numeric(individ_45_adult_ses_t$parameter),as.numeric(neigh_45_adult_ses_t$parameter)),
  p_value    = c(individ_ses_t$p.value,neigh_ses_t$p.value,
                 individ_45_child_ses_t$p.value,neigh_45_child_ses_t$p.value,
                 individ_45_adult_ses_t$p.value,neigh_45_adult_ses_t$p.value),
  conf_low   = c(individ_ses_t$conf.int[1],neigh_ses_t$conf.int[1],
                 individ_45_child_ses_t$conf.int[1],neigh_45_child_ses_t$conf.int[1],
                 individ_45_adult_ses_t$conf.int[1],neigh_45_adult_ses_t$conf.int[1]),
  conf_high  = c(individ_ses_t$conf.int[2],neigh_ses_t$conf.int[2],
                 individ_45_child_ses_t$conf.int[2],neigh_45_child_ses_t$conf.int[2],
                 individ_45_adult_ses_t$conf.int[2],neigh_45_adult_ses_t$conf.int[2])
)

full_t_test_results <- results_list_t_separated_full %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) 

write.csv(full_t_test_results, file = "full_t_test_results.csv", row.names = TRUE)

#-------------------------------------------------------------------------------
# between-group statistical significance between covariates/no covariates (create loop through all combos)

results_list_plot_2 <- list()

counter <- 1

for (t in unique(performance_combined$time_assessed)) {
  for (c in unique(performance_combined$ses_var_rename)) {
    
    subset_data <- subset(performance_combined, 
                          time_assessed == t & ses_var_rename == c)
    
    test <- t.test(r ~ covariates_yn, data = subset_data)
    
    results_list_plot_2[[counter]] <- data.frame(
      time_assessed = t,
      ses_var_rename  = c,
      t_value    = as.numeric(test$statistic),
      df         = as.numeric(test$parameter),
      p_value    = test$p.value,
      conf_low   = test$conf.int[1],
      conf_high  = test$conf.int[2]
    )
    
    counter <- counter + 1
  }
}

results_df_plot_2 <- do.call(rbind, results_list_plot_2)

results_df_covariates <- results_df_plot_2 %>%
  group_by(ses_var_rename) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  ungroup()
write.csv(results_df_covariates, file = "results_df_covariates.csv", row.names = TRUE)

#-------------------------------------------------------------------------------
# Ensure that categorical variables are factors, define levels for plotting
performance_combined$time_assessed <- factor(performance_combined$time_assessed, levels = c("Childhood", "Adulthood", "Age 45"))
performance_combined$time_assessed_rename <- factor(performance_combined$time_assessed_rename, levels = c("Childhood (birth-15)", "Adulthood (26-45)", "Age 45 (concurrent with scan)"))
performance_combined$time_assessed_split <- factor(performance_combined$time_assessed_split, levels = c("Childhood (birth-15)", "Childhood (3-11)","Adulthood (26-45)", "Age 45 (concurrent with scan)"))
performance_combined$ses_var_rename <- factor(performance_combined$ses_var_rename, levels = c("Individual SES","Neighborhood SES"))

#-------------------------------------------------------------------------------
# FDR-adjust p-values from null distribution for multiple comparisons

p_values_nocovar <- c(.148, .229, .021, .318, .045, .218)

p_adjusted_nocovar <- p.adjust(p_values_nocovar, method = "fdr")

p_values_covar <- c(.187, .304, .038, .318, .550, .129,.414)

p_adjusted_covar <- p.adjust(p_values_nocovar, method = "fdr")

#-------------------------------------------------------------------------------
# Visualizations
#-------------------------------------------------------------------------------
df_check_lims <- performance_combined %>% 
  filter(covariates_yn == "N" & time_assessed_rename !="Age 45 (concurrent with scan)")
df_check_lims2 <- performance_combined %>% 
  filter(covariates_yn == "N")
df_check_lims3 <- performance_combined %>% 
  filter(covariates_yn == "N" & time_assessed_rename =="Age 45 (concurrent with scan)")

range(df_check_lims$r, na.rm = T)
range(df_check_lims2$r, na.rm = T)
range(df_check_lims3$r, na.rm = T)
#-------------------------------------------------------------------------------
figure_2 <- ggplot(data = performance_combined %>% 
                               filter(covariates_yn == "N" & time_assessed_rename !="Age 45 (concurrent with scan)" & ses_var_rename == "Individual SES"), 
                             aes(x=time_assessed_split, y=r, fill=time_assessed_split)) + 
  geom_rain(alpha = .5,
            rain.side = "f1x1",
            point.args = list(aes(color = time_assessed_split, alpha = 0.5)))+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#F0E442","#0072B2"
  )) +
  scale_color_manual(values = c("#F0E442", "#0072B2")) +
  ylim(-0.25, 0.5) +
  guides(fill = 'none', color = 'none', alpha = "none") +
  labs(title = "Individual SES and age 45 functional connectivity",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Time point assessed")+
  theme(text = element_text(size = 14)) 
ggsave("figure_2.pdf", plot = figure_2, dpi = 300, width = 7, height = 6, units = "in") 

#-------------------------------------------------------------------------------
figure_S2_A <- ggplot(data = performance_combined %>% 
                             filter(covariates_yn == "N" & ses_var_rename == "Individual SES"), 
                           aes(x=time_assessed_split, y=r, fill=time_assessed_split)) + 
  geom_rain(alpha = .5,
            rain.side = "l",
            point.args = list(aes(color = time_assessed_split, alpha = 0.5)))+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#F0E442","#0072B2","grey"
  )) +
  scale_color_manual(values = c("#F0E442", "#0072B2","grey")) +
  ylim(-0.25, 0.5) +
  guides(fill = 'none', color = 'none', alpha = "none") +
  labs(title = "Prediction performance:\n Individual SES and age 45 functional connectivity",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Time point assessed")+
  theme(text = element_text(size = 14)) 
ggsave("figure_S2_A.pdf", plot = figure_S2_A, dpi = 300, width = 7, height = 6, units = "in") 

#-------------------------------------------------------------------------------
figure_4 <- ggplot(data = performance_combined %>% 
                             filter(covariates_yn == "N" & time_assessed_rename !="Age 45 (concurrent with scan)" & ses_var_rename == "Neighborhood SES"), 
                           aes(x=time_assessed_split, y=r, fill=time_assessed_split)) + 
  geom_rain(alpha = .5,
            rain.side = "f1x1",
            point.args = list(aes(color = time_assessed_split, alpha = 0.5)))+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#F0E442","#0072B2"
  )) +
  scale_color_manual(values = c("#F0E442", "#0072B2")) +
  ylim(-0.25, 0.5) +
  guides(fill = 'none', color = 'none', alpha = "none") +
  labs(title = "Prediction performance:\n Neighborhood SES and age 45 functional connectivity",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Time point assessed")+
  theme(text = element_text(size = 14)) 
ggsave("figure_4.pdf", plot = figure_4, dpi = 300, width = 7, height = 6, units = "in") 

#-------------------------------------------------------------------------------
figure_S2_B <- ggplot(data = performance_combined %>% 
                           filter(covariates_yn == "N" & ses_var_rename == "Neighborhood SES"), 
                         aes(x=time_assessed_split, y=r, fill=time_assessed_split)) + 
  geom_rain(alpha = .5,
            rain.side = "l",
            point.args = list(aes(color = time_assessed_split, alpha = 0.5)))+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#F0E442","#0072B2","grey"
  )) +
  scale_color_manual(values = c("#F0E442", "#0072B2","grey")) +
  ylim(-0.25, 0.5) +
  guides(fill = 'none', color = 'none', alpha = "none") +
  labs(title = "Prediction performance:\n Neighborhood SES and age 45 functional connectivity",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Time point assessed")+
  theme(text = element_text(size = 14)) 
ggsave("figure_S2_B.pdf", plot = figure_S2_B, dpi = 300, width = 7, height = 6, units = "in") 

#-------------------------------------------------------------------------------
covariates_combined_perf <- ggplot(data = performance_combined %>% filter(time_assessed_rename !="Age 45 (concurrent with scan)"), 
                                   aes(x=time_assessed_rename, y=r, fill=covariates_yn)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#af8dc3","#7fbf7b"
  )) +
  labs(title = "Prediction accuracy of socioeconomic variables for age 45 functional connectivity",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Covariates added y/n") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 14)) +
  facet_wrap(~ ses_var_rename)
ggsave("combined_perf_cov_facet.pdf", plot = covariates_combined_perf, dpi = 300, width = 10, height = 7, units = "in") 


performance_combined_cov_var <- performance_combined %>% 
  mutate(ses_detail_name = case_when(ses_var_rename == "Individual SES" ~ "Individual SES, covarying for neighborhood SES (green)",
                                     ses_var_rename == "Neighborhood SES" ~ "Neighborhood SES, covarying for individual SES (green)"))
#-------------------------------------------------------------------------------
covariates_combined_perf_individ <- ggplot(data = performance_combined %>% filter(time_assessed_rename !="Age 45 (concurrent with scan)" & ses_var_rename == "Individual SES"), 
                                           aes(x=time_assessed_rename, y=r, fill=covariates_yn)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(-0.34, 0.5)) +
  scale_fill_manual(values = c("#af8dc3","#7fbf7b"
  )) +
  labs(title = "Prediction accuracy of individual SES \nwithout (purple) and with (green) covarying for neighborhood SES",
       y = "Prediction accuracy (r)", x = " ",
       fill = "Covariates added y/n") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 14),
        legend.position = "none") 

covariates_combined_perf_neigh <- ggplot(data = performance_combined %>% filter(time_assessed_rename !="Age 45 (concurrent with scan)" & ses_var_rename == "Neighborhood SES"), 
                                         aes(x=time_assessed_split, y=r, fill=covariates_yn)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(-0.34, 0.5)) +
  scale_fill_manual(values = c("#af8dc3","#7fbf7b"
  )) +
  labs(title = "Prediction accuracy of neighborhood SES \nwithout (purple) and with (green) covarying for individual SES",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Covariates added y/n") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 14),
        legend.position = "bottom") 

figure_S4 <- ggarrange(covariates_combined_perf_individ, covariates_combined_perf_neigh,
                       ncol = 1,
                       labels = c("A", "B"))

ggsave("figure_S4.pdf", plot = figure_S4, dpi = 300, width = 8, height = 8, units = "in") 
#-------------------------------------------------------------------------------
figure_S5 <- ggplot(data = performance_combined_cov_var %>% filter(time_assessed_rename =="Age 45 (concurrent with scan)"), 
                                   aes(x=time_assessed_rename, y=r, fill=covariates_yn)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#af8dc3","#7fbf7b"
  )) +
  labs(title = "Prediction performance: age 45 socioeconomic status and \nconcurrent (age 45) functional connectivity, without (purple) and with (green) covariates",
       y = "Prediction accuracy (r)", x = "Time of socioeconomic status assessment",
       fill = "Covariates added y/n") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 14),
        legend.position = "bottom") +
  facet_wrap(~ ses_detail_name)
ggsave("figure_S5.pdf", plot = figure_S5, dpi = 300, width = 10, height = 7, units = "in") 


