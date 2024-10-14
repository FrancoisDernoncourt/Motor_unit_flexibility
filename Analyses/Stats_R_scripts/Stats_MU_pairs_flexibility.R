# Loading many packages
# install.packages("olsrr")
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

# Get package citation
citation() #for r
citation("lme4") #for linear model
citation("emmeans") #for linear model

### Initializing parameters

flexibility_metric_to_use <- "max_displacement_normalized_trials_mean" #"cross_corr_trials_mean" #"max_dispersion_normalized_trials_mean" # "max_displacement_normalized_trials_mean"
# Name of the metric to use for flexibility (according to the column names in the .csv generated in Matlab)
# The normalized values are initially present only for the "dispersion" metric,
# but I re-create the other normalized values later in the R script by dividing them by the normalization coefficient specific to the MU pair and condition.
# Normalization coefficient = sum of max DR of both MUs from the pair (nb of norms used for dispersion calculation)
muscle_group <- c("Quadriceps","Triceps_surae")
root_folder <- c() # Root folders to look for the subfolders containing the flexibility table (.csv) of each subject for each condition
root_folder[1] <- "D:/THESE/_DONNEES_MANIP/VL_GM_Compare/Population_flexibility_with_all_quadriceps/QUADRI_flexibility" # for VL
root_folder[2] <- "D:/THESE/_DONNEES_MANIP/VL_GM_Compare/Population_flexibility_with_all_quadriceps/GM_flexibility" # for GM
output_folder <- "D:/THESE/_DONNEES_MANIP/VL_GM_Compare/Population_flexibility_with_all_quadriceps/R_Output" # Output folder
subject_idx_max <- c() # The subjects should be orderly indexed It won't work with subject c(1,2,4) for example, because "3" is missing
subject_idx_max[1] <- 7 # number of subjects for VL
subject_idx_max[2] <- 6 # number of subjects for GM
condition_names <- c("Plateau","0.25hz sinusoids","1hz sinusoids","3hz sinusoids",
                     "all sinusoid conditions","all conditions")
only_MUs_used_for_FA = FALSE # Use only MUs that were included in the factor analysis (MUs which were continuous on the plateau)
only_MUs_matched_in_all_conditions = FALSE # Use only MUs that are present/matched in all conditions. With a few exceptions: some subjects had 0 matched (or even decomposed) MUs in sin3 for the GM, so for these subjects I considered all MUs to be matched if they were present in conditions "plateau", "sin0.25" and "sin1", but not necessarily "sin3".
conditions_to_remove_for_stats = c("Plateau",
                                   "all sinusoid conditions","all conditions") # Remove all these conditions for the stats
muscle_pairs_to_keep <- c("VL-VL","GM-GM")
remove_subjects <- c() #c("50") for example
report_only_mean_of_all_conditions = FALSE
log_transform_values = FALSE # Log-transform the flexibility values or not (this makes their distribution more normal-looking, though it's far from perfect)
# random_effects = c("Subject","MU_pair","None") # All the random effects to try out
random_effects = c("Subject")

# Adding some info about the initializing parameters in the names of the output files so I don't get lost about which file corresponds to which analysis
if (only_MUs_used_for_FA == TRUE)
{
  savefile_suffix_info = "_only_FA_MUs"
} else {
  savefile_suffix_info = ""
}
if (only_MUs_matched_in_all_conditions == TRUE)
{
  savefile_suffix_info = paste(savefile_suffix_info, "_only_MUs_matched_in_all_cond", sep="")
}
if (log_transform_values == TRUE)
{
  savefile_suffix_info = paste(savefile_suffix_info, "_log_transform", sep="")
}

######################
### CREATE THE DATASET
######################

df <- data.frame()
# Load each .csv file
for (muscle_group_i in 1:length(muscle_group))
{
  generic_folder <- paste("_",muscle_group[muscle_group_i],"_MU_pairs",sep="")
  for (subjecti in 1:subject_idx_max[muscle_group_i])
  {
    new_df <- read.csv(paste(root_folder[muscle_group_i],"/S",subjecti, generic_folder,"/",
                             "flexibility_general_table_output.csv", sep=""),
                       header = TRUE, sep = ",")
    new_df$Subject <- subjecti * 10^(muscle_group_i-1) # Assign a specific subject index to each participant
    new_df$Muscle_group <- muscle_group[muscle_group_i]
    
    # if muscle_group == triceps_surae, fill every missing column with GM
    if (muscle_group[muscle_group_i] == "Triceps_surae")
    {
      new_df$MU_x_muscle <- "GM"
      new_df$MU_y_muscle <- "GM"
      new_df$MU_pair_muscles <- "GM-GM"
    }
    
    # Check if MU Y (2nd MU of the pair) was used in FA
    # This should have been saved automatically in the .csv Matlab output, but I made a mistake.
    # In the .csv file, "MU_y_FA_status" = "MU_x_FA_status"
    # So, I am correcting the mistake in the R script directly
    # I am fixing the mistake by re-assigning the correct "MU_y_FA_status". I do so by finding the "MU_x_FA_status" of the same MU.
    new_df$MU_y_FA_status <- -2
    included_in_FA_index_temp <- data.frame(index = new_df$MU_x_matched, FA_status = new_df$MU_x_FA_status)
    included_in_FA_index_temp <- included_in_FA_index_temp[!duplicated(included_in_FA_index_temp), ]
    # Replace "MU_y_FA_status" according to "MU_x_FA_status" of the same corresponding MU
    for (rowi in 1:nrow(new_df))
    {
      MU_y_temp = new_df$MU_y_matched[rowi]
      MU_index_temp <- which(included_in_FA_index_temp$index == MU_y_temp)
      if (length(MU_index_temp) >= 1)
      {
        if (included_in_FA_index_temp$FA_status[MU_index_temp] > 0)
        {
          new_df$MU_y_FA_status[rowi] = included_in_FA_index_temp$FA_status[MU_index_temp]
        }
      }
    }
    # Bind the temporary data frame to the main data frame
    df <- rbind(df,new_df)
  }
}
remove(new_df)

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
df$Muscle_group <- factor(df$Muscle_group)
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

# Get values and display them
df$flexibility <- df[[flexibility_metric_to_use]]
if (log_transform_values == TRUE) # Log-transformed value, if the option has been selected
{
  find_min_value_above_zero <- df %>%
    dplyr::filter(flexibility > 0)
  find_min_value_above_zero <- find_min_value_above_zero$flexibility
  find_min_value_above_zero <- min(find_min_value_above_zero)
  df$flexibility <- log(df$flexibility+find_min_value_above_zero)
}

# Assigning condition names instead of just "1" for platau, "2" for sin 0.25, etc. ...
df$Condition <- factor(df$Condition, levels = levels(droplevels(df$Condition)) )
levels(df$Condition) <- condition_names

# Save data frame as csv
write.csv(df, "flexibility_output.csv", row.names = FALSE)


#############
# DESCRIPTIVE PLOTS OF FLEXIBILITY ACCORDING TO MUSCLE-PAIRING
#############

# Re-writing it here for quick iteration
# flexibility_metric_to_use <- "cross_corr_trials_mean"
# df$flexibility <- df[[flexibility_metric_to_use]]

multiplot_temp <- c()
plot_iter <- 0
scale_max <- max(df$flexibility, na.rm = TRUE)
scale_min <- min(df$flexibility, na.rm = TRUE)
for (conditioni in 1:(length(levels(df$Condition))-2)) # -2 because discarding "all conditions", since very very few RF MUs
{
    plot_iter <- plot_iter+1

    df_to_plot_temp <- df %>%
      filter(Condition == condition_names[conditioni])
    
    multiplot_temp[[plot_iter]] <- ggplot(data = df_to_plot_temp,
                                          aes(x = MU_pair_muscles, y = flexibility, fill = MU_pair_muscles)) +
      geom_boxplot(width = 0.5, size = 1) +
      coord_cartesian(ylim=c(scale_min,scale_max)) +
      stat_summary(aes(group = MU_pair_muscles), # Add mean values to box plots
                   fun = mean, geom = "point", shape = 20, size = 5) +
      labs(title = condition_names[conditioni],
           y = flexibility_metric_to_use,
           x = "") +
      theme_pubr() + #theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))
    
    multiplot_temp[[plot_iter]] <- multiplot_temp[[plot_iter]] +
      scale_x_discrete(drop = FALSE)
}
muscle_pairs_plot <-arrangeGrob(
  grobs = multiplot_temp,
  ncol = length(multiplot_temp))
# Save the "normality check" plot
title_compare_muscle_pairs <- textGrob(
  paste(flexibility_metric_to_use, sep=""),
  gp = gpar(fontsize = 16, fontface = "bold"))
muscle_pairs_plot <- arrangeGrob(muscle_pairs_plot, top = title_compare_muscle_pairs)
grid.draw(muscle_pairs_plot)
ggsave(filename = paste(output_folder,"/",flexibility_metric_to_use,
                        "_according_to_muscles",
                        savefile_suffix_info,".png",
                        sep=""),
       plot = muscle_pairs_plot,
       limitsize = FALSE,
       width = 30,
       height = 10,
       dpi = 100)









###########
### SUBSAMPLE DATASET FOR THE STATISTICAL ANALYSIS
### Not every condition, not every MU pair...
###########

df_to_use <- df %>% filter(MU_pair_muscles %in% muscle_pairs_to_keep)
df_to_use$MU_pair_muscles <- droplevels(df_to_use$MU_pair_muscles)
df_to_use$Muscle <- df_to_use$MU_pair_muscles

if (only_MUs_matched_in_all_conditions) # Selecting only MU pairs which have been matched in all conditions
{
  df_subgroup_temp <- df_to_use[df_to_use$Condition == "all conditions",]
  mu_pairs_idx_matched_in_all_cond <- df_subgroup_temp$MU_pair
  df_to_use <- subset(df_to_use, MU_pair %in% mu_pairs_idx_matched_in_all_cond)
}

if (only_MUs_used_for_FA) # Selecting only MU pairs for which both MUs have been used for the FA on the plateau ( = MUs which were continuous on the plateau)
{
  df_to_use <- subset(df_to_use, MU_x_FA_status == 1)
  df_to_use <- subset(df_to_use, MU_y_FA_status == 1)
  # Code:
  # 1 if the MU was present on the plateau condition, it has been matched at least once, and it was continuous on the plateau (used for FA)
  # 0 if the MU was present on the plateau condition, it has been matched at least once, but it wasn't continuous on the plateau (not used for FA)
  # -1 if the MU was present on the plateau condition, it was continuous (used for FA), but it hasn't been matched in any other condition
  # -2 if the MU wasn't continuous (not used for FA) and it wasn't matched in any other condition (whether present on plateau or present only on another condition)
}

if (length(remove_subjects) >= 1)
{
  for (subjecti in 1:length(remove_subjects))
  {
    df_to_use <- df_to_use[df_to_use$Subject != remove_subjects[subjecti],]
  }
}

df_to_use$flexibility <- df_to_use[[flexibility_metric_to_use]]
# find empty rows becaue of unrepresented conditions
rows_to_remove <- which(is.na(df_to_use$flexibility))
df_to_use <- df_to_use[-rows_to_remove, ]

########
### PLOT FLEXIBILITY ACCORDING TO MUSCLE AND CONDITION (FOR MU PAIRS)
########

signif_symbol <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                      symbols = c("**** p < 0.0001", "*** p < 0.001", "** p < 0.01", "* p < 0.05", "p > 0.05 - ns"))
# ^ This was just to directly put significant differences in the plots, but I am not using it

# Create plot
plot_temp <- ggplot(data = df_to_use,
                    aes(x = Condition, y = flexibility, fill = Muscle)) +
  geom_boxplot(width = 0.5, size = 1) +
  stat_summary(aes(group = interaction(Condition, Muscle), color = Muscle), # Add mean values to box plots
               fun = mean, geom = "point", shape = 20, size = 5,
               position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("firebrick4","royalblue4")) +
  labs(title = paste("Flexibility (",flexibility_metric_to_use,") according to muscles - ",
                     savefile_suffix_info,sep=""),
       y = flexibility_metric_to_use) +
  theme_pubr() + #theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))
# Save plot
ggsave(
  filename = paste(output_folder,"/",
                   flexibility_metric_to_use, savefile_suffix_info,".png",
                   sep=""),
  plot = plot_temp,
  width = 14,
  height = 6,
  dpi = 300)

########
### PLOT FLEXIBILITY ACCORDING TO CONDITION FOR EACH SUBJECT (for MU pairs)
########
multiplot_temp <- c()
plot_iter <- 0
colors_fill <- c("turquoise3","salmon")
colors_mean <- c("royalblue4","firebrick4")
scale_max <- max(df_to_use$flexibility, na.rm = TRUE)
scale_min <- min(df_to_use$flexibility, na.rm = TRUE)
for (muscli in 1:length(muscle_group))
{
  for (subjecti in 1:subject_idx_max[muscli])
  {
    plot_iter <- plot_iter+1
    current_subject_idx_to_plot <- subjecti * 10^(muscli-1)
    # current_subject_idx_to_plot <- factor(current_subject_idx_to_plot)
    df_to_plot_temp <- df_to_use %>%
      filter(Subject == current_subject_idx_to_plot)
    if (nrow(df_to_plot_temp) == 0) { # if plot is empty because subject has been removed
      multiplot_temp[[plot_iter]] <- ggplot()
      next
    }
    
    levels(df_to_plot_temp$Condition) <-
      c("plateau","sin0.25","sin1","sin3","all_sin","all_cond")
    
    multiplot_temp[[plot_iter]] <- ggplot(data = df_to_plot_temp,
      aes(x = Condition, y = flexibility)) +
      geom_boxplot(width = 0.5, size = 1, fill = colors_fill[muscli]) +
      coord_cartesian(ylim=c(scale_min,scale_max)) +
      stat_summary(aes(group = Condition), # Add mean values to box plots
                   fun = mean, geom = "point", shape = 20, size = 5,
                   color = colors_mean[muscli]) +
      labs(title = paste(muscle_group[muscli], " - Subject ",subjecti,sep=""),
           y = flexibility_metric_to_use,
           x = "") +
      theme_pubr() + #theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))
      
      multiplot_temp[[plot_iter]] <- multiplot_temp[[plot_iter]] +
        scale_x_discrete(drop = FALSE)
  }
}
per_subject_plot <-grid.arrange(
  grobs = multiplot_temp,
  ncol = max(subject_idx_max))
# Save the "normality check" plot
ggsave(filename = paste(output_folder,"/",flexibility_metric_to_use,
                        "Individual subjects data",
                        savefile_suffix_info,".png",
                        sep=""),
       plot = per_subject_plot,
       limitsize = FALSE,
       width = 35,
       height = 10,
       dpi = 100)

##############
### STATISTICAL TESTS
##############

output_folder_stat <- paste(output_folder,"/Stats/",sep="")
dir.create(output_folder_stat)

# Removing conditions for statistical tests
for (conditioni in 1:length(conditions_to_remove_for_stats))
{
  df_to_use <- df_to_use[df_to_use$Condition != conditions_to_remove_for_stats[conditioni], ]
}
if (report_only_mean_of_all_conditions == TRUE)
{
  temp_mean_flexibility_over_conditions <- df_to_use %>%
    group_by(MU_pair) %>% summarize(flexibility = mean(flexibility))
  df_to_use$Condition <- "mean for all sinusoid conditions"
  df_to_use$Condition <- as.factor(df_to_use$Condition)
  df_to_use <- df_to_use %>%
    distinct(MU_pair, .keep_all = TRUE)
  df_to_use <- merge(df_to_use,
                     temp_mean_flexibility_over_conditions[, c("MU_pair", "flexibility")], by = "MU_pair")
  # two values are created for flexibility = flexibility.x (from df_to_use) and flexibility.y (from temp_mean_flexibility_over_conditions)
  df_to_use$flexibility = df_to_use$flexibility.y
  df_to_use$flexibility.x <- NULL
  df_to_use$flexibility.y <- NULL
}
df_to_use$Condition <- factor(df_to_use$Condition, levels = levels(droplevels(df_to_use$Condition)) )

# Check normality of the resulting flexibility values (but not necesary for the linear model)
plot_title <- "QQ plot of flexibility values"
if (log_transform_values == TRUE)
{
  plot_title <- paste(plot_title," - Log transformed", sep="")
}
plot_qq <- c()
# Histogram (look at distribution). The "binwidth" has to be changed according to the flexibility value selected
plot_qq[[1]] <- ggplot(df_to_use, aes(x = flexibility)) + 
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = plot_title, x = "Flexibility values", y = "Frequency")
# QQ plot comparing to a normal distribution
plot_qq[[2]] <- ggplot(df_to_use, aes(sample = flexibility)) +
  geom_qq() +
  geom_qq_line(color = "red", linewidth = 1) +
  labs(title = "QQ Plot")
plots_qq <-grid.arrange(
  grobs = plot_qq,
  ncol = 2)
# Save the "normality check" plot
ggsave(filename = paste(output_folder_stat,
                        flexibility_metric_to_use,
                        savefile_suffix_info,"_",
                        "data_distribution plot",".png",
                        sep=""),
       plot = plots_qq,
       width = 10,
       height = 5,
       dpi = 100)

########
### PLOT FLEXIBILITY ACCORDING TO CONDITION FOR EACH SUBJECT (for MU pairs) - ACTUAL DATA USED
### (SO WITH CONDITIONS TO BE REMOVED ACTUALLY REMOVED)
########

colors_fill <- c("turquoise3","salmon")
colors_mean <- c("royalblue4","firebrick4")
multiplot_temp <- c()
multiplot_temp[[1]] <- ggplot(data = df_to_use,
                    aes(x = reorder(Subject, flexibility, mean),
                        y = flexibility, fill = Muscle)) +
  geom_boxplot(width = 0.5, size = 1) +
  stat_summary(aes(group = Subject), # Add mean values to box plots
               fun = mean, geom = "point", shape = 20, size = 5) +
  labs(title = "Flexibility for each subject, independent of condition",
       y = flexibility_metric_to_use,
       x = "Subject") +
  theme_pubr() + #theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))

multiplot_temp[[2]] <- ggplot(data = df_to_use,
                              aes(x = reorder(Muscle, flexibility, mean), y = flexibility, fill = Muscle)) +
  geom_boxplot(width = 0.5, size = 1) +
  stat_summary(aes(group = Muscle), # Add mean values to box plots
               fun = mean, geom = "point", shape = 20, size = 5) +
  labs(title = "Flexibility per muscle, all subjects aggregated",
      y = flexibility_metric_to_use,
       x = "Muscle") +
  theme_pubr() + #theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))

per_subject_plot <- grid.arrange(multiplot_temp[[1]], multiplot_temp[[2]],
                                 ncol = 2, widths = c(3, 1))

ggsave(filename = paste(output_folder,"/",flexibility_metric_to_use,
                        "_data_for_stats",
                        savefile_suffix_info,".png",
                        sep=""),
       plot = per_subject_plot,
       limitsize = TRUE,
       width = 15,
       height = 8,
       dpi = 100)

############
### CREATE THE LINEAR MODELS
############

linear_model_no_random_effect <- lm(flexibility ~ Muscle * Condition,
                                    data = df_to_use)
linear_model_simple_random_effect <- lmer(flexibility ~ Muscle * Condition +
                       (1 | Subject), data = df_to_use)
linear_model_complex_random_effect <- lmer(flexibility ~ Muscle * Condition +
              (Condition | Subject), data = df_to_use)

compare_models <- anova(linear_model_simple_random_effect,
      linear_model_no_random_effect,
      linear_model_complex_random_effect)

write_xlsx(tidy(compare_models),
           path = paste(output_folder_stat,
                        flexibility_metric_to_use,
                        savefile_suffix_info,
                        "_compare_lmms.xls",sep=""))

linear_model <- linear_model_complex_random_effect

# anova() can work with a mixture lmer and lm models. However, due to the way that R's type system is set up, it only works if the first argument is an lmer model. That is, the anova.merMod() method (which gets called if the first argument is a [g]lmer model) knows how to deal with lm objects, but the anova.lm() method (which gets called if an lm object is first) doesn't know about merMod objects ...


anova_results <- anova(linear_model)
# Save output as .xls
write_xlsx(tidy(anova_results),
           path = paste(output_folder_stat,
                        flexibility_metric_to_use,
                        savefile_suffix_info,
                        "_anova.xls",sep=""))


summary(linear_model)
model_predictions <- lmerTest::ls_means(linear_model, confit = 0.95)
model_predictions$Name <- rownames(model_predictions)
write_xlsx(model_predictions,
           path = paste(output_folder_stat,
                        flexibility_metric_to_use,
                        savefile_suffix_info,
                        "_model_predictions.xls",sep=""))
# Fixed effects
plot_to_save <- emmip(linear_model, ~ Muscle, CIs = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_fixed_effect_muscle.png",sep=""),
       plot_to_save, width=6, height=6, dpi=100)
plot_to_save <- emmip(linear_model, ~ Condition, CIs = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_fixed_effect_condition.png",sep=""),
       plot_to_save, width=6, height=6, dpi=100)
plot_to_save <- emmip(linear_model, ~ Muscle:Condition, CIs = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_fixed_effect_interaction.png",sep=""),
       plot_to_save, width=6, height=6, dpi=100)
# Post-hoc stratified fixed effects
plot_to_save <- emmip(linear_model, ~ Muscle | Condition, CIs = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_post_hoc_muscle_by_condition.png",sep=""),
       plot_to_save, width=6, height=6, dpi=100)
plot_to_save <- emmip(linear_model, ~ Condition | Muscle, CIs = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_post_hoc_condition_by_muscle.png",sep=""),
       plot_to_save, width=6, height=6, dpi=100)

## Post-hoc comparisons
emmeans_obj <- emmeans(linear_model, ~ Muscle | Condition)
emmansplot <- plot(emmeans_obj, comparisons = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_post_hoc_compare_muscles.png",sep=""),
       emmansplot, width=6, height=6, dpi=100)
pwc_within_conditions <- pairs(emmeans_obj, adjust = "bonf") # Comparing muscles, stratifying by condition
eff_size(emmeans_obj, sigma = sigma(linear_model), edf = Inf)
# Save output as .xls
write_xlsx(tidy(pwc_within_conditions),
           path = paste(output_folder_stat,
                        flexibility_metric_to_use,
                        savefile_suffix_info,
                        "_post_hoc_compare_muscles.xls",sep=""))

emmeans_obj <- emmeans(linear_model, ~ Condition | Muscle, repeated = "Subject")
emmansplot <- plot(emmeans_obj, comparisons = TRUE)
ggsave(paste(output_folder_stat,
             flexibility_metric_to_use,
             savefile_suffix_info,
             "_post_hoc_compare_conditions.png",sep=""),
       emmansplot, width=6, height=6, dpi=100)
pwc_within_muscle <- pairs(emmeans_obj, adjust = "bonf") # Comparing conditions, stratifying by muscle
# Save output as .xls
write_xlsx(tidy(pwc_within_muscle),
           path = paste(output_folder_stat,
                        flexibility_metric_to_use,
                        savefile_suffix_info,
                        "_post_hoc_compare_conditions.xls",sep=""))


##### PLOT THE VALUES USED FOR NORMALIZAION
# Correspond to the mean (for each trial) of the sum of the max firing rates of both MUs
plot_temp <- ggplot(data = df_to_use,
                    aes(x = Condition, y = normalization_value_trials_mean, fill = Muscle)) +
  geom_boxplot(width = 0.5, size = 1) +
  stat_summary(aes(group = interaction(Condition, Muscle), color = Muscle),
               fun = mean, geom = "point", shape = 20, size = 5,
               position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("firebrick4","royalblue4")) +
  labs(title = paste("Sum of max DRs (mean over trials) according to condition and muscle - ",
                     savefile_suffix_info,sep="")) +
  theme_pubr() + #theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))

ggsave(
  filename = paste(output_folder,"/normalization_values",
                   "_sum_of_max_DRs_", savefile_suffix_info,".png",
                   sep=""),
  plot = plot_temp,
  width = 14,
  height = 6,
  dpi = 300)

###### PERFORM LINEAR MODEL VISUAL CHECKS

# check_model() takes a long time, so I commented out everything

### Only on the main model of interest (linear model on MU pairs, with subject as random effect)
  # linear_model <- lmer(flexibility ~ Muscle + Condition +
  #                        (1 | Subject), data = df_to_use)
# summary(linear_model)
# anova(linear_model)
# ### Check model assumptions and and assess model quality
# ### Visual check of model various assumptions (normality of residuals, normality of random effects, heteroscedasticity, homogeneity of variance, multicollinearity).
plot_model_check <- check_model(linear_model)
png(file = paste(output_folder_stat,
                 flexibility_metric_to_use,
                 savefile_suffix_info,
                 "_model_checks.png",sep=""),
    width=2000, height=2000)
plot_model_check
dev.off()
# ### Checks for and locates influential observations (i.e., "outliers") via several distance and/or clustering methods.
# check_outliers(linear_model)
# ### Posterior predictive checks mean "simulating replicated data under the fitted model and then comparing these to the observed data" 
plot_predictions <- check_predictions(linear_model)
png(file = paste(output_folder_stat,
                 flexibility_metric_to_use,
                 savefile_suffix_info,
                 "_predictions.png",sep=""),
    width=1000, height=700)
plot_predictions
dev.off()

###############
# PLOT RESIDUALS
###############

residual_plot <- c()
resid <- data.frame(residuals(linear_model))
density_temp = density(resid$residuals.linear_model.)
scale_max_y <- ceiling(max(density_temp$y, na.rm = TRUE))

residual_plot[[1]] <- ggplot(resid, aes(x=residuals.linear_model.)) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill="pink") +
  geom_density(alpha = 0.3, fill = "tomato") +
  labs(title = paste("Distribution of residuals ",sep=""),
       y = "density (%)") +
  theme_pubr() + #theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))
hist_resid_data <- ggplot_build(residual_plot[[1]])$data[[1]]
scale_min_x <- min(hist_resid_data$x)
scale_max_x <- max(hist_resid_data$x)

residual_plot[[1]] <- residual_plot[[1]] +
  coord_cartesian(ylim=c(0,scale_max_y),
                  xlim=c(scale_min_x,scale_max_x)) 

resid_mean <- mean(residuals(linear_model))
resid_sd <- sd(residuals(linear_model))
theoritical_normal_distrib <- data.frame(rnorm(nrow(resid), resid_mean, resid_sd))
  
residual_plot[[2]] <- ggplot(theoritical_normal_distrib,
                             aes(x=rnorm.nrow.resid...resid_mean..resid_sd.)) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill="skyblue") +
  geom_density(alpha = 0.3, fill = "dodgerblue") +
  coord_cartesian(ylim=c(0,scale_max_y),
                  xlim=c(scale_min_x,scale_max_x)) +
  labs(title = paste("Theoritical normal distribution ",sep=""),
       y = "density (%)") +
  theme_pubr() + #theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))

residuals_mutliplot <- arrangeGrob(residual_plot[[1]], residual_plot[[2]],
                                 ncol = 2)
kolmogorv_smirnov_test_result <- ks.test(resid$residuals.linear_model., "pnorm")
if (kolmogorv_smirnov_test_result$p.value < 0.05) {
  distrib_text = " (the residuals are unlikely to come from a normal distribution)"
} else {
  distrib_text = " (the distribution of residuals is compatible with what would be expected from a normal distribution)"
}
title_residuals_plot <- textGrob(
  paste("Kolmogorov-Smirnov test p-value = ",
      kolmogorv_smirnov_test_result$p.value,
      distrib_text,  sep=""),
  gp = gpar(fontsize = 12, fontface = "bold"))
residuals_mutliplot <- arrangeGrob(residuals_mutliplot, top = title_residuals_plot)
grid.draw(residuals_mutliplot)

ggsave(
  filename = paste(output_folder_stat,
                   flexibility_metric_to_use,
                   "_residuals_histogram",savefile_suffix_info
                   ,".png",
                   sep=""),
  plot = residuals_mutliplot,
  width = 10,
  height = 6,
  dpi = 100)

stop("END OF SCRIPT")

# ############### OTHER CHECKS
# 
# # Checking other stuffs
# ggplot(linear_model, aes(.fitted, .resid)) + geom_point() +
#   facet_wrap(~ Muscle + Condition, ncol = 4)
# 
# plot(linear_model, resid(.) ~ fitted(.) | Muscle + Condition )
# 
# linear_model <- lmer(log(flexibility) ~ Muscle * Condition +
#                        (1 | Subject), data = df_to_use)
# 
# # Counting MUs for plateau
# test <- df %>% filter(Condition=="Plateau") %>% group_by(Subject) %>%
#   distinct(MU_x) %>% count(Subject)
# mean(test$n[c(1,2,3,4,5,6,7)]) # VL
# sd(test$n[c(1,2,3,4,5,6,7)]) # VL
# mean(test$n[c(8,9,10,11,12,13)]) # GM
# sd(test$n[c(8,9,10,11,12,13)]) # GM

# SOME DESCRIPTIVE VALUES FOR THE PAPER

df_only_ref_condition <- df %>%
  filter(Condition=="3hz sinusoids") #(Condition=="Plateau")
  # filter(Condition!="all sinusoid conditions" &
  #          Condition!="all conditions") #(Condition=="Plateau")
for (muscle_groupi in 0:6)
{
  if (muscle_groupi == 0) {
    temp_output <- df_only_ref_condition %>%
      filter(MU_pair_muscles=="RF-VL" |
               MU_pair_muscles=="RF-VM")
    print(paste("RF-VL & RF-VM - Mean dispersion = ",round(mean(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
                " ; Std Dev = ",round(sd(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
                sep=""))
    print(paste("RF-VL & RF-VM - Mean displacement = ",round(mean(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
                " ; Std Dev = ",round(sd(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
                sep=""))
    next
  } else if (muscle_groupi == 1) {
    MU_pair_muscle_to_use <- "RF-VL"
  } else if (muscle_groupi == 2) {
    MU_pair_muscle_to_use <- "RF-VM"
  } else if (muscle_groupi == 3) {
    MU_pair_muscle_to_use <- "VL-VL"
  } else if (muscle_groupi == 4) {
    MU_pair_muscle_to_use <- "GM-GM"
  } else if (muscle_groupi == 5) {
    MU_pair_muscle_to_use <- "VL-VM"
  } else if (muscle_groupi == 6) {
    MU_pair_muscle_to_use <- "VM-VM"
  }
temp_output <- df_only_ref_condition %>%
  filter(MU_pair_muscles==MU_pair_muscle_to_use)
print(paste(MU_pair_muscle_to_use," - Mean dispersion = ",round(mean(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
            " ; Std Dev = ",round(sd(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
            sep=""))
print(paste(MU_pair_muscle_to_use," - Mean displacement = ",round(mean(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
            " ; Std Dev = ",round(sd(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
            sep=""))
}

# flexibility values for the sinusoidal conditions
for (conditioni in 2:(length(condition_names)-2) )
{
  temp_output <- df_to_use %>%
    filter(Condition==condition_names[conditioni])
  print(paste(condition_names[conditioni]," - Mean dispersion = ",round(mean(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),3),
              " ; Std Dev = ",round(sd(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),3),
              sep=""))
  print(paste(condition_names[conditioni]," - Mean displacement = ",round(mean(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),3),
              " ; Std Dev = ",round(sd(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),3),
              sep=""))
}


########## VALUES ACROSS MUSCLE AND CONDITIONS
for (conditioni in 1:(length(condition_names)-2) )
{
  print(paste(condition_names[conditioni]," : ",sep=""))
  for (muscli in 1:2)
  {
    muscle_pair_to_pick <- c()
    if (muscli == 1) {
      muscle_pair_to_pick <- 'VL-VL'
    } else if (muscli == 2) {
      muscle_pair_to_pick <- 'GM-GM'
    }
    temp_output <- df %>%
      filter(Condition==condition_names[conditioni]) %>%
      filter(MU_pair_muscles==muscle_pair_to_pick)
    print(paste("     ",muscle_pair_to_pick," : ",sep=""))
    print(paste("      - ",condition_names[conditioni]," - Mean dispersion = ",round(mean(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
                  " ; Std Dev = ",round(sd(temp_output$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
                  sep=""))
    print(paste("      - ",condition_names[conditioni]," - Mean displacement = ",round(mean(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
                " ; Std Dev = ",round(sd(temp_output$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
                sep=""))
    print(paste("      - ",condition_names[conditioni]," - Mean correl = ",round(mean(temp_output$cross_corr_trials_mean, na.rm = TRUE),2),
                " ; Std Dev = ",round(sd(temp_output$cross_corr_trials_mean, na.rm = TRUE),2),
                sep=""))
  }
}


#########
# STD of participant's flexibility
#########

df_compare_subject_variability <- data.frame(matrix(ncol = 14, nrow = 0))
colnames(df_compare_subject_variability)<- c('Subject','Muscle',
                                             'sin 0.25 dispersion mean','sin 0.25 dispersion std dev',
                                             'sin 0.25 displacement mean','sin 0.25 displacement std dev',
                                             'sin 1 dispersion mean','sin 1 dispersion std dev',
                                             'sin 1 displacement mean','sin 1 displacement std dev',
                                             'sin 3 dispersion mean','sin 3 dispersion std dev',
                                             'sin 3 displacement mean','sin 3 displacement std dev')
row_iter <- 0
for (subjecti in 1:length(levels(temp_df$Subject)))
{
  row_iter <- row_iter +1
  muscle_pair_to_pick <- c()
  if (subjecti <= 7) {
    muscle_pair_to_pick <- 'VL-VL'
  } else {
    muscle_pair_to_pick <- 'GM-GM'
  }
  
  temp_df <- df_to_use %>%
    filter(Subject==levels(temp_df$Subject)[subjecti]) %>%
    filter(MU_pair_muscles==muscle_pair_to_pick)
  
  print(paste("Subject #", levels(temp_df$Subject)[subjecti], " (",muscle_pair_to_pick,") : ", sep = ""))
  
  df_compare_subject_variability[row_iter,] <- NA
  df_compare_subject_variability$Subject[row_iter] = levels(temp_df$Subject)[subjecti]
  df_compare_subject_variability$Muscle[row_iter] = muscle_pair_to_pick
    
    for (conditioni in 2:(length(condition_names)-2) )
    {
      temp_df_condition <- temp_df %>%
        filter(Condition==condition_names[conditioni])
      
      print(paste("     ",condition_names[conditioni]," :",sep=""))
      
      print(paste("       - dispersion: ", round(mean(temp_df_condition$max_dispersion_normalized_trials_mean, na.rm = TRUE),2),
                  " ± ", round(sd(temp_df_condition$max_dispersion_normalized_trials_mean, na.rm = TRUE),2), sep = ""))
      
      print(paste("       - displacement: ", round(mean(temp_df_condition$max_displacement_normalized_trials_mean, na.rm = TRUE),2),
                  " ± ", round(sd(temp_df_condition$max_displacement_normalized_trials_mean, na.rm = TRUE),2), sep = ""))
      
      condition_first_col_idx <- ((conditioni-2)*4)+3
      df_compare_subject_variability[row_iter,condition_first_col_idx] <- mean(temp_df_condition$max_dispersion_normalized_trials_mean, na.rm = TRUE)
      df_compare_subject_variability[row_iter,condition_first_col_idx+1] <- sd(temp_df_condition$max_dispersion_normalized_trials_mean, na.rm = TRUE)
      df_compare_subject_variability[row_iter,condition_first_col_idx+2] <- mean(temp_df_condition$max_displacement_normalized_trials_mean, na.rm = TRUE)
      df_compare_subject_variability[row_iter,condition_first_col_idx+3] <- sd(temp_df_condition$max_displacement_normalized_trials_mean, na.rm = TRUE)
      
    }
}

df_compare_subject_variability_VL <- df_compare_subject_variability %>%
  filter(Muscle=='VL-VL')
df_compare_subject_variability_GM <- df_compare_subject_variability %>%
  filter(Muscle=='GM-GM')

#       # VL
print(paste(
  'Std dev of mean dispersion of participants for the VL, sin 0.25 condition = ',
  round(sd(df_compare_subject_variability_VL$`sin 0.25 dispersion mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean displacement of participants for the VL, sin 0.25 condition = ',
  round(sd(df_compare_subject_variability_VL$`sin 0.25 displacement mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean dispersion of participants for the VL, sin 1 condition = ',
  round(sd(df_compare_subject_variability_VL$`sin 1 dispersion mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean displacement of participants for the VL, sin 1 condition = ',
  round(sd(df_compare_subject_variability_VL$`sin 1 displacement mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean dispersion of participants for the VL, sin 3 condition = ',
  round(sd(df_compare_subject_variability_VL$`sin 3 dispersion mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean displacement of participants for the VL, sin 3 condition = ',
  round(sd(df_compare_subject_variability_VL$`sin 3 displacement mean`, na.rm = TRUE),2),
  sep = ""))

#       # GM
print(paste(
  'Std dev of mean dispersion of participants for the GM, sin 0.25 condition = ',
  round(sd(df_compare_subject_variability_GM$`sin 0.25 dispersion mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean displacement of participants for the GM, sin 0.25 condition = ',
  round(sd(df_compare_subject_variability_GM$`sin 0.25 displacement mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean dispersion of participants for the GM, sin 1 condition = ',
  round(sd(df_compare_subject_variability_GM$`sin 1 dispersion mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean displacement of participants for the GM, sin 1 condition = ',
  round(sd(df_compare_subject_variability_GM$`sin 1 displacement mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean dispersion of participants for the GM, sin 3 condition = ',
  round(sd(df_compare_subject_variability_GM$`sin 3 dispersion mean`, na.rm = TRUE),2),
  sep = ""))
print(paste(
  'Std dev of mean displacement of participants for the GM, sin 3 condition = ',
  round(sd(df_compare_subject_variability_GM$`sin 3 displacement mean`, na.rm = TRUE),2),
  sep = ""))





test <- data.frame(matrix(ncol = 0, nrow = 100))
test$one <- rnorm(100)
test$two <- rnorm(100)+1
sd( c( mean(test$one) , mean(test$two) ) )
mean( c( sd(test$one) , sd(test$two) ) )
