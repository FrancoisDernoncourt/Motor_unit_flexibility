# Loading many packages
# install.packages("scales")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggnewscale)
library(writexl)
library(gridExtra)
library(grid)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(performance)
library(lmerTest)
library(Matrix)
library(see)
library(patchwork)
library(ggpmisc)
library(broom)
library(MASS)
library(olsrr)
library(scales)

root_folder_flexibility <- "D:/THESE/_DONNEES_MANIP/VL_GM_Compare/Population_flexibility_with_all_quadriceps/QUADRI_flexibility/"
foldername_flexibility_pet_subjet <- "_Quadriceps_MU_pairs/"
filename_to_load_flexibility <- "flexibility_general_table_output.csv"
root_folder_distance <- "D:/THESE/_DONNEES_MANIP/VL/VL_Distance_for_each_MU_pair/"
filename_to_load_distance <- "_VL_distances_between_MUs.csv"
subject_idx_max <- 7
condition_names <- c("Plateau","0.25hz sinusoids","1hz sinusoids","3hz sinusoids",
                     "all sinusoid conditions","all conditions") # consider only conditions 1-4

new_folder <- c( "D:/THESE/_DONNEES_MANIP/VL/VL_Distance_for_each_MU_pair/R_output")
dir.create(new_folder)

# CREATE DATA FRAME
df <- data.frame()
for (subjecti in 1:subject_idx_max) {
  # Load flexibility values
  df_temp_flexibility <- read.csv(paste(root_folder_flexibility,"S",subjecti, foldername_flexibility_pet_subjet, filename_to_load_flexibility, sep=""),
                                  header = TRUE, sep = ",")
  df_temp_flexibility <- df_temp_flexibility %>% filter(Condition < 5)
  df_temp_flexibility <- df_temp_flexibility %>% filter(MU_pair_muscles=="VL-VL")
  # Load distance values
  df_temp_distances <- read.csv(paste(root_folder_distance,"S",subjecti, filename_to_load_distance, sep=""),
                                  header = TRUE, sep = ",")
  # Create temp df from them
  df_temp <- df_temp_flexibility
  df_temp$Distance <- df_temp_distances$Distance
  # Bind the temporary data frame to the main data frame
  df <- rbind(df,df_temp)
}
remove(df_temp_distances, df_temp_flexibility, df_temp)

# Assign index of MU pairs and convert "NaNs" values (from Matlab) into "NA" values (for R)
for (rowi in 1:nrow(df))
{
  #df$MU_pair[rowi] <- paste(df$MU_x_matched[rowi],"-",df$MU_y_matched[rowi],sep="")
  if (df$MU_x_matched[rowi] == "NaN") {
    df$MU_pair[rowi] <- "NaN"
  } else {
    if (df$MU_x_matched[rowi] != 0 && df$MU_y_matched[rowi] != 0)
    {
      # Assign an individual index to each specific pair of MUs
      temp_highest_MU_idx <- max(c(df$MU_x_matched[rowi], df$MU_y_matched[rowi]))
      temp_lowest_MU_idx <- min(c(df$MU_x_matched[rowi], df$MU_y_matched[rowi]))
      df$MU_pair[rowi] <- df$Subject[rowi]*1e4 + temp_lowest_MU_idx*1e2 + temp_highest_MU_idx
    } else {
      df$MU_pair[rowi] <- "NaN"
    }
  }
  for (coli in 1:ncol(df))
  {
    # Change all NaN values to NA
    if (df[rowi,coli] == "NaN")
    {
      df[rowi,coli] = NA
    }
  }
}
# Convert variables into factors
df$Subject <- factor(df$Subject)
df$Condition <- factor(df$Condition)
df$MU_x_muscle <- factor(df$MU_x_muscle)
df$MU_y_muscle <- factor(df$MU_y_muscle)
df$MU_pair_muscles <- factor(df$MU_pair_muscles)
df$MU_pair <- factor(df$MU_pair)

# Get normalization value for the other flexibility metric
# ...for the trials mean values
df$normalization_value_trials_mean <- df$mean_dispersion_trials_mean / df$mean_dispersion_normalized_trials_mean
df$mean_displacement_normalized_trials_mean <- df$mean_displacement_trials_mean / df$normalization_value_trials_mean
df$max_displacement_normalized_trials_mean <- df$max_displacement_trials_mean / df$normalization_value_trials_mean
df$max_dispersion_normalized_trials_mean <- df$max_dispersion_trials_mean / df$normalization_value_trials_mean
# ...for the trials concatenated values
df$normalization_value_trials_concatenated <- df$mean_dispersion / df$mean_dispersion_normalized
df$mean_displacement_normalized <- df$mean_displacement / df$normalization_value_trials_concatenated
df$max_displacement_normalized <- df$max_displacement / df$normalization_value_trials_concatenated
df$max_dispersion_normalized <- df$max_dispersion / df$normalization_value_trials_concatenated

# Set subject colors #need 7 colors
cscale <- c("coral","coral1","tomato","tomato2",
            "indianred2","indianred3","firebrick4")

# Cross-corr according to distance ###################

df_for_plot <- df %>% filter(Condition==1)
# df_for_plot$Subject <- as.numeric(df_for_plot$Subject)
distance_VS_xcorr_plot <- ggplot(data = df_for_plot,
                                         aes(x = Distance, y = cross_corr_trials_mean,
                                             color = Subject)) +
  geom_point2(size = 3.5, stroke = 0, alpha = 0.1) +
  labs(title = "Distance of MUs VS cross-correlation of MUs (Plateau condition)",
       y = "Cross-correl",
       x = "Distance (mm)") +
  geom_line(stat = "smooth", method = lm,
            aes(x = Distance, y = cross_corr_trials_mean, color = Subject),
            size = 1.5, alpha = 0.5) +
  geom_line(stat = "smooth", method = lm,
            aes(x = Distance, y = cross_corr_trials_mean, color=c("pooled")),
            color = 'firebrick2', size = 2, alpha = 1) +
  scale_color_manual(values=cscale) +
  theme_classic2() + 
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
distance_VS_xcorr_plot
ggsave(filename = paste(new_folder,"/distance_VS_xcorr_plot.png",sep=""),
        plot = distance_VS_xcorr_plot)

distance_results <- df_for_plot %>%
  group_by(Subject) %>%
  summarize(
    correl_xcorr_per_subject = cor(Distance, cross_corr_trials_mean),
    rsquared_xcorr_per_subject = summary(lm(cross_corr_trials_mean ~ Distance))$r.squared,
    pvalue_xcorr_per_subject = summary(lm(cross_corr_trials_mean ~ Distance))$coefficients[2,4]
  )
distance_results_all_subjects <- df %>%
  summarize(
    correl_xcorr_per_subject = cor(Distance, cross_corr_trials_mean),
    rsquared_xcorr_per_subject = summary(lm(cross_corr_trials_mean ~ Distance))$r.squared,
    pvalue_xcorr_per_subject = summary(lm(cross_corr_trials_mean ~ Distance))$coefficients[2,4]
  )
distance_results_all_subjects$pvalue_xcorr_per_subject


distance_results$Subject <- as.character(distance_results$Subject)
distance_results[subject_idx_max+1,] = 
  data.frame(Subject = "pooled",
  correl_xcorr_per_subject = cor(df_for_plot$Distance, df_for_plot$cross_corr_trials_mean),
  rsquared_xcorr_per_subject = summary(lm(cross_corr_trials_mean ~ Distance, data = df_for_plot))$r.squared)
distance_results$Subject <- as.factor(distance_results$Subject)

# Dispersion according to distance ###################
multiplot_flexibility <- c()
scale_max_displacement <- max(df$max_displacement_normalized_trials_mean, na.rm = TRUE)
scale_min_displacement <- min(df$max_displacement_normalized_trials_mean, na.rm = TRUE)
scale_max_dispersion <- max(df$max_dispersion_normalized_trials_mean, na.rm = TRUE)
scale_min_dispersion <- min(df$max_dispersion_normalized_trials_mean, na.rm = TRUE)
for (conditioni in 2:4) {
  df_for_plot <- df %>% filter(Condition==conditioni)
  # Dispersion
  multiplot_flexibility[[conditioni-1]] <- ggplot(data = df_for_plot,
                                             aes(x = Distance, y = max_dispersion_normalized_trials_mean,
                                                 color = Subject)) +
      geom_point2(size = 3.5, stroke = 0, alpha = 0.1) +
      labs(title = paste(condition_names[conditioni]," - Dispersion",sep=""),
           y = "Dispersion",
           x = "Distance (mm)") +
      geom_line(stat = "smooth", method = lm,
                aes(x = Distance, y = max_dispersion_normalized_trials_mean, color = Subject),
                size = 1.5, alpha = 0.5) +
      geom_line(stat = "smooth", method = lm,
                aes(x = Distance, y = max_dispersion_normalized_trials_mean, color=c("pooled")),
                color = 'firebrick2', size = 2, alpha = 1) +
      coord_cartesian(ylim=c(scale_min_dispersion,scale_max_dispersion)) +
      scale_color_manual(values=cscale) +
      theme_classic2() + 
      theme(axis.line = element_line(linewidth = 1, color = "darkgrey"),
            legend.position = "none")
  
  distance_results_temp <- df_for_plot %>%
    group_by(Subject) %>%
    summarize(
      correl_dispersion_per_subject = cor(Distance, max_dispersion_normalized_trials_mean),
      rsquared_dispersion_per_subject = summary(lm(max_dispersion_normalized_trials_mean ~ Distance))$r.squared,
      pval_dispersion_per_subject =  summary(lm(max_dispersion_normalized_trials_mean ~ Distance))$coefficients[2,4]
    )
  distance_results_temp$Subject <- as.character(distance_results_temp$Subject)
  distance_results_temp[subject_idx_max+1,] = 
    data.frame(Subject = "pooled",
               correl_dispersion_per_subject = cor(df_for_plot$Distance, df_for_plot$max_dispersion_normalized_trials_mean),
               rsquared_dispersion_per_subject = summary(lm(max_dispersion_normalized_trials_mean ~ Distance, data = df_for_plot))$r.squared,
               pval_dispersion_per_subject =  summary(lm(max_dispersion_normalized_trials_mean ~ Distance, data = df_for_plot))$coefficients[2,4]
               )
  distance_results_temp$Subject <- as.factor(distance_results_temp$Subject)
  
  distance_results <- cbind(distance_results, distance_results_temp)
  
  # Displacement
  multiplot_flexibility[[3+(conditioni-1)]]<- ggplot(data = df_for_plot,
                                                aes(x = Distance, y = max_displacement_normalized_trials_mean,
                                                    color = Subject)) +
      geom_point2(size = 3.5, stroke = 0, alpha = 0.1) +
      labs(title = paste(condition_names[conditioni]," - Displacement",sep=""),
           y = "Displacement",
           x = "Distance (mm)") +
      geom_line(stat = "smooth", method = lm,
                aes(x = Distance, y = max_displacement_normalized_trials_mean, color = Subject),
                size = 1.5, alpha = 0.5) +
      geom_line(stat = "smooth", method = lm,
                aes(x = Distance, y = max_displacement_normalized_trials_mean, color=c("pooled")),
                color = 'firebrick2', size = 2, alpha = 1) +
      coord_cartesian(ylim=c(scale_min_displacement,scale_max_displacement)) +
      scale_color_manual(values=cscale) +
      theme_classic2() + 
      theme(axis.line = element_line(linewidth = 1, color = "darkgrey"),
            legend.position = "none")
  
  distance_results_temp <- df_for_plot %>%
    group_by(Subject) %>%
    summarize(
      correl_displacement_per_subject = cor(Distance, max_displacement_normalized_trials_mean),
      rsquared_displacement_per_subject = summary(lm(max_displacement_normalized_trials_mean ~ Distance))$r.squared,
      pval_displacement_per_subject =  summary(lm(max_displacement_normalized_trials_mean ~ Distance))$coefficients[2,4]
      )
  distance_results_temp$Subject <- as.character(distance_results_temp$Subject)
  distance_results_temp[subject_idx_max+1,] = 
    data.frame(Subject = "pooled",
               correl_displacement_per_subject = cor(df_for_plot$Distance, df_for_plot$max_displacement_normalized_trials_mean),
               rsquared_displacement_per_subject = summary(lm(max_displacement_normalized_trials_mean ~ Distance, data = df_for_plot))$r.squared,
               pval_displacement_per_subject =  summary(lm(max_displacement_normalized_trials_mean ~ Distance, data = df_for_plot))$coefficients[2,4]
               )
  distance_results_temp$Subject <- as.factor(distance_results_temp$Subject)
  
  distance_results <- cbind(distance_results, distance_results_temp)
}
distance_VS_flexibility_plot <-arrangeGrob(
  grobs = multiplot_flexibility,
  ncol = 3)
grid.draw(distance_VS_flexibility_plot)
ggsave(filename = paste(new_folder,"/distance_VS_flexibility_plots.png",sep=""),
       plot = distance_VS_flexibility_plot, limitsize = FALSE, width = 20, height = 10, dpi = 100)

write_xlsx(distance_results, path = paste(new_folder,"/distance_results.xls",sep=""))

# distance_results <- subset(distance_results, select = -c("Subject.1", "Subject.2", "Subject.3", "Subject.4", "Subject.5", "Subject.6"))
# distance_results <- distance_results %>% rename(correl_dispersion_per_subject = correl_dispersion_per_subject_sin0_25,
#                                                 rsquared_dispersion_per_subject = rsquared_dispersion_per_subject_sin0_25,
#                                                 correl_displacement_per_subject = correl_displacement_per_subject_sin0_25,
#                                                 rsquared_displacement_per_subject = rsquared_displacement_per_subject_sin0_25,
#                                                 
#                                                 correl_dispersion_per_subject.1 = correl_dispersion_per_subject_sin1,
#                                                 rsquared_dispersion_per_subject.1 = rsquared_dispersion_per_subject_sin1,
#                                                 correl_displacement_per_subject.1 = correl_displacement_per_subject_sin1,
#                                                 rsquared_displacement_per_subject.1 = rsquared_displacement_per_subject_sin1,
#                                                 
#                                                 correl_dispersion_per_subject.2 = correl_dispersion_per_subject_sin3,
#                                                 rsquared_dispersion_per_subject.2 = rsquared_dispersion_per_subject_sin3,
#                                                 correl_displacement_per_subject.2 = correl_displacement_per_subject_sin3,
#                                                 rsquared_displacement_per_subject.2 = rsquared_displacement_per_subject_sin3)


###### statistical tests
stat_model <- lmer(max_displacement_normalized_trials_mean ~ Distance + (Condition | Subject), data = df)
summary(stat_model)
fixed_effects <- fixef(stat_model)
fixed_effects_pvalues <- coef(summary(stat_model))
fixed_effects_df <- as.data.frame(fixed_effects_pvalues)

stat_model <- lmer(max_dispersion_normalized_trials_mean ~ Distance + (Condition | Subject), data = df)
summary(stat_model)
fixed_effects <- fixef(stat_model)
fixed_effects_pvalues <- coef(summary(stat_model))
fixed_effects_df <- as.data.frame(fixed_effects_pvalues)

