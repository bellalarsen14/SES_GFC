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
library(forcats)

################### OPTIONS #############################
root <- "/parentfoldername/"
setwd("filepathforyourfigures")

# define the behavioral (SES) variables for analysis
behavvar_list <- data.frame(behavvar = c("ADI311","SESchildhd",
                                         "PH45_AreaDeptot","SESall45",
                                         "neighdep2645_factor","ses_composite"))

# load performance statistics
performance_combined <- read.csv(paste0(root,'performance_combined_update.csv')) 

# load behavioral dataframe
behavdata <- read.csv(file = paste0(root,"person_level_df_all_ses_comp.csv"))

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
# between-age statistical significance (create loop through all combos)
individ_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% filter(covariates_yn == "N" & time_assessed != "Age 45" & ses_var_rename == "Individual SES"))

neigh_ses_t <- t.test(r ~ time_assessed, data = performance_combined %>% filter(covariates_yn == "N" & time_assessed != "Age 45" & ses_var_rename == "Neighborhood SES"))

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
# create a long version of behavioral data, merge with functional connectivity
behavdata_merge <- behavdata %>% 
  select(snum,sex, neighdep2645_factor,ses_composite,SESchildhd,ADI311,SESall45,PH45_AreaDeptot, AverageFD)

gfc_snums <- GFC %>% 
  pull(snum)

# filter for only people with functional connectivity data - the study sample
behavdata_merged_fc <- behavdata_merge %>% 
  filter(snum %in% gfc_snums) %>% 
  select(snum, ADI311, SESchildhd, PH45_AreaDeptot, SESall45, neighdep2645_factor, ses_composite, sex, AverageFD)

# create z-scored versions of the SES scores for visualization, invert neighborhood SES
behavdata_merged_fc_scaled <- behavdata_merged_fc %>%
  mutate(ADI311_flip = -1*ADI311,
         PH45_AreaDeptot_flip = -1*PH45_AreaDeptot,
         neighdep2645_factor_flip = -1*neighdep2645_factor) %>% 
  mutate(ADI311_z = as.numeric(scale(ADI311_flip)),
         SESchildhd_z = as.numeric(scale(SESchildhd)),
         PH45_AreaDeptot_z = as.numeric(scale(PH45_AreaDeptot_flip)),
         SESall45_z = as.numeric(scale(SESall45)),
         neighdep2645_factor_z = as.numeric(scale(neighdep2645_factor_flip)),
         ses_composite_z = as.numeric(scale(ses_composite)))

# Sample characteristics
behavdata_merged_fc_scaled %>% 
  count(sex) %>% 
  mutate(perc = n/sum(n))

vars <- c("ADI311_flip", "neighdep2645_factor_flip", "PH45_AreaDeptot_flip",
          "PH45_AreaDeptot_z", "SESchildhd", "ses_composite", "SESall45")

sapply(behavdata_merged_fc_scaled[vars], function(x) c(mean = mean(x, na.rm = TRUE),
                                                       sd = sd(x, na.rm = TRUE)))

# correlate childhood and adulthood individual SES
cor.test(behavdata_merged_fc_scaled$SESchildhd_z, behavdata_merged_fc_scaled$ses_composite_z, method = "pearson")

# correlate childhood and adulthood neighborhood SES
cor.test(behavdata_merged_fc_scaled$ADI311_z, behavdata_merged_fc_scaled$neighdep2645_factor_z, method = "pearson")

# create a long version of behavioral data for visualization
long_behavdata_merged <- pivot_longer(
  data = behavdata_merged_fc_scaled,
  cols = c(SESchildhd_z, ADI311_z, ses_composite_z, neighdep2645_factor_z,
           SESall45_z,PH45_AreaDeptot_z),   # columns to gather
  names_to = "ses_var",        # new column for old column names
  values_to = "z_score_value"           # new column for values
) %>% 
  mutate(measure_level = case_when(ses_var == "ADI311_z" | ses_var == "neighdep2645_factor_z" | ses_var == "PH45_AreaDeptot_z" ~ "Neighborhood SES",
                                   ses_var == "SESchildhd_z" | ses_var == "ses_composite_z" | ses_var == "SESall45_z" ~ "Individual SES")) %>% 
  mutate(time_assessed = case_when(ses_var == "ADI311_z" | ses_var == "SESchildhd_z" ~ "Childhood (birth-15)",
                                   ses_var == "neighdep2645_factor_z" | ses_var == "ses_composite_z" ~ "Adulthood (26-45)",
                                   ses_var == "PH45_AreaDeptot_z" | ses_var == "SESall45_z" ~ "Age 45 (concurrent w/ scan)"))

long_behavdata_merged$ses_var <- factor(long_behavdata_merged$ses_var, levels = c("SESchildhd_z", "ADI311_z", 
                                                                                  "ses_composite_z","neighdep2645_factor_z",
                                                                                  "SESall45_z", "PH45_AreaDeptot_z"))
long_behavdata_merged$measure_level <- factor(long_behavdata_merged$measure_level, levels = c("Individual SES","Neighborhood SES"))
long_behavdata_merged$time_assessed <- factor(long_behavdata_merged$time_assessed, levels = c("Childhood (birth-15)","Adulthood (26-45)","Age 45 (concurrent w/ scan)"))

long_behavdata_merged <- long_behavdata_merged %>% 
  mutate(time_assessed_noage = case_when(time_assessed == "Childhood (birth-15)" ~ "Childhood",
                                         time_assessed == "Adulthood (26-45)" ~ "Adulthood",
                                         time_assessed == "Age 45 (concurrent w/ scan)" ~ "Age 45 (concurrent w/ scan)")) 

long_behavdata_merged$time_assessed_noage <- factor(long_behavdata_merged$time_assessed_noage, levels = c("Childhood", "Adulthood", "Age 45 (concurrent w/ scan)"))

#-------------------------------------------------------------------------------
# Visualizations: main text
#-------------------------------------------------------------------------------
# Check the range of values for standardized axes

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
# Figure 2

figure_1 <- long_behavdata_merged %>%
  filter(time_assessed != "Age 45 (concurrent w/ scan)") %>% 
  ggplot( aes(x=z_score_value, color=time_assessed_noage, fill=time_assessed_noage)) +
  geom_histogram(alpha=0.6, binwidth = 1) +
  scale_fill_manual(values = c("#F0E442","#0072B2"
  )) +
  scale_color_manual(values = c("#F0E442","#0072B2"
  )) +
  theme_classic() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14)
  ) +
  facet_grid(vars(measure_level),vars(time_assessed_noage))+
  labs(title = "Distribution of socioeconomic status measures in childhood and adulthood",
       x = "Socioeconomic status (z-score)", y = "Count", fill = "SES measure", color = "SES measure")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 16))

ggsave("figure_1.pdf", plot = figure_1, dpi = 300, width = 10, height = 6, units = "in") 

#-------------------------------------------------------------------------------
# Figure 2

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
# note: Figure 3 is in an additional R file, along with feature importance analyses
#-------------------------------------------------------------------------------

# Figure 4

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
# Visualizations: supplement 
#-------------------------------------------------------------------------------
# Figure S1
figure_S1 <- long_behavdata_merged %>% 
  ggplot( aes(x=z_score_value, color=time_assessed_noage, fill=time_assessed_noage)) +
  geom_histogram(alpha=0.6, binwidth = 1) +
  scale_fill_manual(values = c("#F0E442","#0072B2","grey"
  )) +
  scale_color_manual(values = c("#F0E442","#0072B2","grey"
  )) +
  theme_classic() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14)
  ) +
  facet_grid(vars(measure_level),vars(time_assessed_noage))+
  labs(title = "Distribution of socioeconomic status measures in childhood and adulthood",
       x = "Socioeconomic status (z-score)", y = "Count", fill = "SES measure", color = "SES measure")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 16))

ggsave("figure_S1.pdf", plot = figure_S1, dpi = 300, width = 12, height = 6, units = "in") 

#-------------------------------------------------------------------------------

# Figure S2
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
#note: Figure S3 is in an additional R file, along with feature importance analyses
#-------------------------------------------------------------------------------

# Figure S4

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
# Figure S5
  # create a detailed name for plotting
performance_combined_cov_var <- performance_combined %>% 
  mutate(ses_detail_name = case_when(ses_var_rename == "Individual SES" ~ "Individual SES, covarying for neighborhood SES (green)",
                                     ses_var_rename == "Neighborhood SES" ~ "Neighborhood SES, covarying for individual SES (green)"))

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

#-------------------------------------------------------------------------------
# note: Figures S6 and S7 are in an additional R file, along with feature importance analyses
