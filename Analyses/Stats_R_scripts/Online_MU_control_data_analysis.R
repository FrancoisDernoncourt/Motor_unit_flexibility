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

nb_subjects <- 6
muscles <- c("VL","GM")
folder <- "D:/THESE/_DONNEES_MANIP/Online/To_use_for_R/"
filename_suffix <- "_Table_output.csv"
min_RT_diff <- 3 #2

df <- data.frame()
for (subjecti in 1:nb_subjects) {
  for (muscli in 1:length(muscles)) {
    current_muscle <- muscles[muscli]
    path_temp <- paste(folder,"S",subjecti,"_",current_muscle,filename_suffix,sep="")
    if (file.exists(path_temp)) {
      df_temp <- read.csv(path_temp,
               header = TRUE, sep = ",")
    } else {
      next
    }
    df <- rbind(df, df_temp)
  }
}
remove(df_temp)

# Convert "NaNs" values (from Matlab) into "NA" values (for R)
for (rowi in 1:nrow(df))
{
  for (coli in 1:ncol(df))
  {
    # Change all NaN values to NA
    if (!is.na(df[rowi,coli])) {
      if (df[rowi,coli] == "NaN")
      {
        df[rowi,coli] = NA
      }
    } 
  }
}

df$Subject <- as.factor(df$Subject)
df$Muscle <- as.factor(df$Muscle)
df$MUx_position <- as.factor(df$MUx_position)
df$MUy_position <- as.factor(df$MUy_position)
df$Target <- as.factor(df$Target)
df$Target_violating_size_principle <- as.factor(df$Target_violating_size_principle)
df$condition <- df$Target=='baseline'
df$condition <- factor(df$condition, levels = c(TRUE,FALSE),
                      labels = c("Baseline", "Dissociation_attempt"))
df$MUpair_unique_idx <- paste(df$Subject,"-",df$Muscle,"_",df$MU_pair)
df$MUpair_unique_idx <- as.factor(df$MUpair_unique_idx)
df$Time_percent_success_rate[df$Displayed_as_feedback==0] <- NA
# Get unique idx of MU pairs displayed as feedback
MUs_index_displayed_as_feedback <- 
  df$MUpair_unique_idx[df$Displayed_as_feedback==1]
MUs_index_displayed_as_feedback <- droplevels(MUs_index_displayed_as_feedback)
MUs_index_displayed_as_feedback <- unique(MUs_index_displayed_as_feedback)

df$Current_target_imply_size_principle_violation <-
  as.character(df$Target) == as.character(df$Target_violating_size_principle)
df$Current_target_imply_size_principle_violation <- as.factor(
  df$Current_target_imply_size_principle_violation)

###############################
## QUANTITATIVE DESCRIPTIONS OF MUs
###############################

df_subset <- df %>% filter(Displayed_as_feedback == 1)

nb_displayed_MUs_VL_per_subject <- df_subset %>% filter(Muscle=="VL") %>% count(Subject, .drop = TRUE)
nb_displayed_MUs_VL_per_subject <- (nb_displayed_MUs_VL_per_subject$n)/2
total_MU_pairs_displayed_VL <- sum(nb_displayed_MUs_VL_per_subject)
mean(nb_displayed_MUs_VL_per_subject)
sd(nb_displayed_MUs_VL_per_subject)

nb_displayed_MUs_GM_per_subject <- df_subset %>% filter(Muscle=="GM") %>% count(Subject, .drop = TRUE)
nb_displayed_MUs_GM_per_subject <- (nb_displayed_MUs_GM_per_subject$n)/2
total_MU_pairs_displayed_GM <- sum(nb_displayed_MUs_GM_per_subject)
mean(nb_displayed_MUs_GM_per_subject)
sd(nb_displayed_MUs_GM_per_subject)

# count all MUs
MU_per_participant_VL <- df %>% filter(Muscle=="VL") %>% group_by(Subject) %>% count(MUx, .drop = TRUE) %>%
  count(Subject, .drop = TRUE)
MU_per_participant_VL <- MU_per_participant_VL$n + 1 # MUx, so need +1 to account for the last MU (which is in MUy)
mean(MU_per_participant_VL)
sd(MU_per_participant_VL)

MU_per_participant_GM <- df %>% filter(Muscle=="GM") %>% group_by(Subject) %>% count(MUx, .drop = TRUE) %>%
  count(Subject, .drop = TRUE)
MU_per_participant_GM <- MU_per_participant_GM$n + 1 # MUx, so need +1 to account for the last MU (which is in MUy)
mean(MU_per_participant_GM)
sd(MU_per_participant_GM)

################################
# VIOLATION OF SIZE PRINCIPLE ==
################################

df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff)
# df_baseline_vs_control <- aggregate(
#   Time_percent_violation ~ 
#     (condition * Subject * Muscle * MUpair_unique_idx),
#   data = df_subset, FUN = mean, na.rm = FALSE)

df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Muscle == "VL") %>%
  filter(condition=="Baseline")
mean(df_subset$Time_percent_violation, na.rm = TRUE)
sd(df_subset$Time_percent_violation, na.rm = TRUE)
min(df_subset$Time_percent_violation, na.rm = TRUE)
max(df_subset$Time_percent_violation, na.rm = TRUE)
df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Muscle == "VL") %>%
  filter(condition!="Baseline")
mean(df_subset$Time_percent_violation, na.rm = TRUE)
sd(df_subset$Time_percent_violation, na.rm = TRUE)
min(df_subset$Time_percent_violation, na.rm = TRUE)
max(df_subset$Time_percent_violation, na.rm = TRUE)

df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Muscle == "GM") %>%
  filter(condition=="Baseline")
mean(df_subset$Time_percent_violation, na.rm = TRUE)
sd(df_subset$Time_percent_violation, na.rm = TRUE)
min(df_subset$Time_percent_violation, na.rm = TRUE)
max(df_subset$Time_percent_violation, na.rm = TRUE)
df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Muscle == "GM") %>%
  filter(condition!="Baseline")
mean(df_subset$Time_percent_violation, na.rm = TRUE)
sd(df_subset$Time_percent_violation, na.rm = TRUE)
min(df_subset$Time_percent_violation, na.rm = TRUE)
max(df_subset$Time_percent_violation, na.rm = TRUE)

### APPROPRIATE MODEL
lm_violation_all_pairs.model <- lmer(
  Time_percent_violation ~ condition + (Muscle * condition | Subject),
  data =  df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff))
summary(lm_violation_all_pairs.model)
lm_violation_all_pairs.p_val <- anova(lm_violation_all_pairs.model)

### Just t-test to check difference between proportion for GM and VL
df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff) %>%
  filter(condition!="Baseline")
t.test(Time_percent_violation ~ Muscle,
       data = df_subset, paired = FALSE)
# Check histogram o eyeball data normality
# removing data = 0 because overwhelmingly numerous
df_subset <- df_subset %>% filter(Time_percent_violation > 0.1)
df_subset$logval = log(df_subset$Time_percent_violation)
ggplot(df_subset, aes(x = logval, fill = Muscle)) +
  geom_histogram(color = "white", alpha = 0.5) +
  labs(title = "Histogram", x = "% time violation", y = "Count") +
  # xlim(c(-1,25)) +
  theme_minimal()

wilcox.test(logval ~ Muscle,
            data = df_subset, paired = FALSE)


# STAT MODELS
check_models.no_random <- lm(Time_percent_violation ~ Muscle * condition,
                      data =  df_subset)
check_models.subject_intercept <- lmer(Time_percent_violation ~ Muscle * condition +
                            (1 | Subject), data =  df_subset)
check_models.subject_cond_on_vars <- lmer(Time_percent_violation ~ Muscle * condition +
                                       (Muscle * condition | Subject), data =  df_subset)
check_models.random_but_no_muscle <- lmer(Time_percent_violation ~ condition +
                                                (condition | Subject), data =  df_subset)
check_models.subject_cond_on_vars_RT_diff <- lmer(Time_percent_violation ~ Muscle * condition * RT_diff_betwwen_MUs +
                                            (Muscle * condition | Subject), data =  df_subset) 
check_models.subject_cond_on_vars_RT_diff_no_muscle <- lmer(Time_percent_violation ~ condition * RT_diff_betwwen_MUs +
                                                    (condition | Subject), data =  df_subset) 
# testing out stuff, most of it is not relevant
anova(check_models.subject_intercept, 
      check_models.subject_cond_on_vars,
      check_models.no_random,
      check_models.random_but_no_muscle) #,
      # check_models.subject_cond_on_vars_RT_diff,
      # check_models.subject_cond_on_vars_RT_diff_no_muscle)

summary(check_models.subject_cond_on_vars_RT_diff_no_muscle)

mean(df_subset$Time_percent_violation) # Mean % time violating size principle
sd(df_subset$Time_percent_violation)
df_temp <- df_subset %>% filter(condition=='Baseline')
mean(df_temp$Time_percent_violation)
sd(df_temp$Time_percent_violation)
df_temp <- df_subset %>% filter(condition=='Dissociation_attempt')
mean(df_temp$Time_percent_violation, na.rm = TRUE)
sd(df_temp$Time_percent_violation, na.rm = TRUE)

# Simple model (but not lowest AIC and BIC values)
test_model <-lmer(Time_percent_violation ~ condition +
                    (condition | Subject), data =  df_subset)
summary(test_model)
# ^ significant with (condition | Subject) & extremely significant with (1 | Subject)
# Talk about it in a comment in the paper, maybe...
# Funny stuff can happen when "wiggling" (or RT evaluation inaccuracy), but not volitionally controlable (success rate is 0)
# Limits = screen for MUs which are most likely to be "dissociable" (Formento)
# (right now was based on filter quality, so more of a pragmatic criterion)
# size principle reversal was rare (give mean), but it did happen on some MUs (up to 100% in some trials). No significiant difference despite heterogeneity, but suggests that flexibility is available, just not volitionally controllable
# I don't talk about it but I also have all the data about grid distal or proximal and target upper (distal) or lower (proximal), maybe something to look into also

lm_violation_all_pairs.model <- lmer(
  Time_percent_violation ~ condition + (Muscle * condition | Subject),
  data =  df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff))
summary(lm_violation_all_pairs.model)
lm_violation_all_pairs.p_val <- anova(lm_violation_all_pairs.model)
# # #
emmeans_obj <- emmeans(lm_violation_all_pairs.model, ~ Muscle)
pwc_muscles_violation_size_principle <- pairs(emmeans_obj, adjust = "bonf")
lm_violation_all_pairs.effsize_muscle <- eff_size(emmeans_obj, sigma =
                                        sigma(lm_violation_all_pairs.model),
                                      edf = df.residual(lm_violation_all_pairs.model))
# # #
emmeans_obj <- emmeans(lm_violation_all_pairs.model, ~ condition)
pwc_condition_violation_size_principle <- pairs(emmeans_obj, adjust = "bonf")
lm_violation_all_pairs.effsize_condition <- eff_size(emmeans_obj, sigma =
                                             sigma(lm_violation_all_pairs.model),
                                           edf = df.residual(lm_violation_all_pairs.model))

plot_size_principle_violation <- ggplot(data = df_subset,
                                        aes(x = condition, y = Time_percent_violation,
                                            color = Muscle, fill = Muscle)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               alpha = 0.4, linewidth = 1.5, width = 0.5) +
  geom_jitter2(position = position_jitterdodge(
    dodge.width = 0.6, jitter.width = 0.53),
    size = 4, alpha = 0.25) +
  scale_color_manual(values=c("dodgerblue","tomato")) +
  scale_fill_manual(values=c("dodgerblue","tomato")) +
  labs(x = "Condition", y = "% Time violation of size principle", 
       title = paste("Size principle violation ; all MU pairs
       Only MU pairs with RT diff > ",min_RT_diff,"% MVC",sep=""), font = "Roboto") +
  theme_classic2() + 
  ylim(c(0,100)) +
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
plot_size_principle_violation
ggsave(paste(folder,"Size_principle_violation_all_MUs.png",sep=""))

###############################################################
# VIOLATION OF SIZE PRINCIPLE == Only MUs displayed as feedback
###############################################################

df_subset_only_displayed_MUs <- df_baseline_vs_control %>%
  filter(MUpair_unique_idx %in% MUs_index_displayed_as_feedback)
# convert MUpair idx to characters and then factors (I don't know why but the factors' levels went to 1)
df_subset_only_displayed_MUs$MUpair_unique_idx <- as.character(
  df_subset_only_displayed_MUs$MUpair_unique_idx)
df_subset_only_displayed_MUs$MUpair_unique_idx <- as.factor(
  df_subset_only_displayed_MUs$MUpair_unique_idx)

# STAT MODELS
lm_violation_only_feedback_pairs.model <- lmer(
  Time_percent_violation ~ Muscle * condition + (Muscle*condition | Subject),
  data =  df_subset_only_displayed_MUs)
summary(lm_violation_only_feedback_pairs.model)
lm_violation_only_feedback_pairs.p_val <- anova(lm_violation_only_feedback_pairs.model)
# # #
emmeans_obj <- emmeans(lm_violation_only_feedback_pairs.model, ~ Muscle)
pwc_muscles_violation_size_principle <- pairs(emmeans_obj, adjust = "bonf")
lm_violation_only_feedback_pairs.effsize_muscle <- eff_size(emmeans_obj, sigma =
                                   sigma(lm_violation_only_feedback_pairs.model),
                                   edf = df.residual(lm_violation_only_feedback_pairs.model))
# # #
emmeans_obj <- emmeans(lm_violation_only_feedback_pairs.model, ~ condition)
pwc_condition_violation_size_principle <- pairs(emmeans_obj, adjust = "bonf")
lm_violation_only_feedback_pairs.effsize_condition <- eff_size(emmeans_obj, sigma =
                                                              sigma(lm_violation_only_feedback_pairs.model),
                                                            edf = df.residual(lm_violation_only_feedback_pairs.model))

plot_size_principle_violation_only_feedback_MUS <- ggplot(
  data = df_subset_only_displayed_MUs,
  aes(x = condition, y = Time_percent_violation,
      color = Muscle, fill = Muscle)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               alpha = 0.4, linewidth = 1.5, width = 0.5) +
  
  geom_jitter(data = df_subset_only_displayed_MUs %>% filter(Muscle=="GM"),
              aes(x = condition, y = Time_percent_violation,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                x = -0.15, nudge.from = c("jittered"), seed = 256),
              size = 5, alpha = 0.6) +
  geom_path(data = df_subset_only_displayed_MUs %>% filter(Muscle=="GM"),
            aes(x = condition, y = Time_percent_violation,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
              x = -0.15, nudge.from = c("jittered"), seed = 256),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  geom_jitter(data = df_subset_only_displayed_MUs %>% filter(Muscle=="VL"),
              aes(x = condition, y = Time_percent_violation,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                x = 0.15, nudge.from = c("jittered"), seed = 256),
              size = 5, alpha = 0.6) +
  geom_path(data = df_subset_only_displayed_MUs %>% filter(Muscle=="VL"),
            aes(x = condition, y = Time_percent_violation,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
              x = 0.15, nudge.from = c("jittered"), seed = 256),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
 
  scale_color_manual(values=c("dodgerblue","tomato")) +
  scale_fill_manual(values=c("dodgerblue","tomato")) +
  labs(x = "Condition", y = "% Time violation of size principle", 
       title = "Only MUs for which feedback was displayed
       Only MU pairs with RT diff > 2% MVC", font = "Roboto") +
  theme_classic2() + 
  ylim(c(0,100)) + 
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
plot_size_principle_violation_only_feedback_MUS
ggsave(paste(folder,"Size_principle_violation_only_feedback_MUs.png",sep=""))

################################
# VIOLATION OF SIZE PRINCIPLE ==
# ONLY MUs DISPLAYED AS FEEDBACK
# ONLY WHEN ATTEMPTING TO DISSOCIATE IN THE DIRECTION WHICH WOULD VIOLAT THE SIZE PRINCIPLE
################################

df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff) %>%
  filter(Current_target_imply_size_principle_violation=='TRUE')
df_baseline_vs_control <- aggregate(
  Time_percent_violation ~ 
    (condition * Subject * Muscle * MUpair_unique_idx),
  data = df_subset, FUN = mean, na.rm = FALSE)

df_subset_only_displayed_MUs <- df_baseline_vs_control %>%
  filter(MUpair_unique_idx %in% MUs_index_displayed_as_feedback)
# convert MUpair idx to characters and then factors (I don't know why but the factors' levels went to 1)
df_subset_only_displayed_MUs$MUpair_unique_idx <- as.character(
  df_subset_only_displayed_MUs$MUpair_unique_idx)
df_subset_only_displayed_MUs$MUpair_unique_idx <- as.factor(
  df_subset_only_displayed_MUs$MUpair_unique_idx)

# STAT MODELS
lm_violation_only_feedback_pairs.model <- lmer(
  Time_percent_violation ~ Muscle * condition + (Muscle*condition | Subject),
  data =  df_subset_only_displayed_MUs)
summary(lm_violation_only_feedback_pairs.model)
lm_violation_only_feedback_pairs.p_val <- anova(lm_violation_only_feedback_pairs.model)
# # #
emmeans_obj <- emmeans(lm_violation_only_feedback_pairs.model, ~ Muscle)
pwc_muscles_violation_size_principle <- pairs(emmeans_obj, adjust = "bonf")
lm_violation_only_feedback_pairs.effsize_muscle <- eff_size(emmeans_obj, sigma =
                                                              sigma(lm_violation_only_feedback_pairs.model),
                                                            edf = df.residual(lm_violation_only_feedback_pairs.model))
# # #
emmeans_obj <- emmeans(lm_violation_only_feedback_pairs.model, ~ condition)
pwc_condition_violation_size_principle <- pairs(emmeans_obj, adjust = "bonf")
lm_violation_only_feedback_pairs.effsize_condition <- eff_size(emmeans_obj, sigma =
                                                                 sigma(lm_violation_only_feedback_pairs.model),
                                                               edf = df.residual(lm_violation_only_feedback_pairs.model))
plot_size_principle_violation_only_feedback_MUS_target_violates_size_principle <-
  ggplot(
  data = df_subset_only_displayed_MUs,
  aes(x = condition, y = Time_percent_violation,
      color = Muscle, fill = Muscle)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               alpha = 0.4, linewidth = 1.5, width = 0.5) +
  
  geom_jitter(data = df_subset_only_displayed_MUs %>% filter(Muscle=="GM"),
              aes(x = condition, y = Time_percent_violation,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                                              x = -0.15, nudge.from = c("jittered"), seed = 256),
              size = 5, alpha = 0.6) +
  geom_path(data = df_subset_only_displayed_MUs %>% filter(Muscle=="GM"),
            aes(x = condition, y = Time_percent_violation,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
                                            x = -0.15, nudge.from = c("jittered"), seed = 256),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  geom_jitter(data = df_subset_only_displayed_MUs %>% filter(Muscle=="VL"),
              aes(x = condition, y = Time_percent_violation,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                                              x = 0.15, nudge.from = c("jittered"), seed = 256),
              size = 5, alpha = 0.6) +
  geom_path(data = df_subset_only_displayed_MUs %>% filter(Muscle=="VL"),
            aes(x = condition, y = Time_percent_violation,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
                                            x = 0.15, nudge.from = c("jittered"), seed = 256),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  scale_color_manual(values=c("dodgerblue","tomato")) +
  scale_fill_manual(values=c("dodgerblue","tomato")) +
  labs(x = "Condition", y = "% Time violation of size principle", 
       title = "Only MUs for which feedback was displayed
       Only when the displayed target involved violating the size principle
       Only MU pairs with RT diff > 2% MVC", font = "Roboto") +
  theme_classic2() + 
  ylim(c(0,100)) + 
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
plot_size_principle_violation_only_feedback_MUS_target_violates_size_principle
ggsave(paste(folder,"Size_principle_violation_only_feedback_MUs_target_violates_size_principle.png",sep=""))

##################################
## SUCCESS RATE ##################
##################################

lm_sucess_rate.model <- lmer(Time_percent_success_rate ~ Muscle + (Muscle | Subject),
                       data =  df %>% filter(Target != 'baseline'))
summary(lm_sucess_rate.model)
lm_sucess_rate.p_val <- anova(lm_sucess_rate.model)
emmeans_obj <- emmeans(lm_sucess_rate.model, ~ Muscle)
pwc_muscles_succes_rate <- pairs(emmeans_obj, adjust = "bonf")
lm_sucess_rate.effsize <- eff_size(emmeans_obj, sigma = sigma(lm_sucess_rate.model),
                               edf = df.residual(lm_sucess_rate.model))

df_only_GM <- df %>% filter(Muscle=="GM")
df_only_VL <- df %>% filter(Muscle=="VL")

plot_success_rate <- ggplot(data = df %>% filter(Target != 'baseline'),
                      aes(x = Target, y = Time_percent_success_rate,
                          color = Muscle, fill = Muscle)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               alpha = 0.4, linewidth = 1.5, width = 0.5) +
  
  geom_jitter(data = df_only_GM %>% filter(Target != 'baseline'),
              aes(x = Target, y = Time_percent_success_rate,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                x = -0.15, nudge.from = c("jittered"), seed = 123),
              size = 5, alpha = 0.6) +
  geom_path(data = df_only_GM %>% filter(Target != 'baseline'),
            aes(x = Target, y = Time_percent_success_rate,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
              x = -0.15, nudge.from = c("jittered"), seed = 123),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  geom_jitter(data = df_only_VL %>% filter(Target != 'baseline'),
              aes(x = Target, y = Time_percent_success_rate,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                x = 0.15, nudge.from = c("jittered"), seed = 123),
              size = 5, alpha = 0.6) +
  geom_path(data = df_only_VL %>% filter(Target != 'baseline'),
            aes(x = Target, y = Time_percent_success_rate,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
              x = 0.15, nudge.from = c("jittered"), seed = 123),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  
  scale_color_manual(values=c("dodgerblue","tomato")) +
  scale_fill_manual(values=c("dodgerblue","tomato")) +
  labs(x = "Target", y = "% Success rate", font = "Roboto",
       title = "Success rate") +
  scale_x_discrete(label = c("Lower quadrant","Upper quadrant")) +
  theme_classic2() + 
  ylim(c(0,100)) +
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
plot_success_rate
ggsave(paste(folder,"Succes_rate_all_trials_both_target_quadrants.png",sep=""))



##################################
## SUCCESS RATE ##################
### WHEN THE TARGET INVOLVES VIOLATING THE SIZE PRINCIPLE OR NOT
##################################

df_subset <- df %>% filter(Displayed_as_feedback == 1) %>% filter(Target != 'baseline') %>% filter(Muscle=="VL")
succes_rate_VL <- df_subset$Time_percent_success_rate
mean(succes_rate_VL)
sd(succes_rate_VL)

min_RT_diff <- 3
df_subset <- df %>% filter(Displayed_as_feedback == 1) %>% filter(Target != 'baseline') %>% filter(Muscle=="VL") %>%
  filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Current_target_imply_size_principle_violation=='FALSE')
succes_rate_VL_compatible_size_principle <- df_subset$Time_percent_success_rate
mean(succes_rate_VL_compatible_size_principle)
sd(succes_rate_VL_compatible_size_principle)

df_subset <- df %>% filter(Displayed_as_feedback == 1) %>% filter(Target != 'baseline') %>% filter(Muscle=="VL") %>%
  filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Current_target_imply_size_principle_violation=='TRUE')
succes_rate_VL_uncompatible_size_principle <- df_subset$Time_percent_success_rate
mean(succes_rate_VL_uncompatible_size_principle)
sd(succes_rate_VL_uncompatible_size_principle)

df_subset <- df %>% filter(Displayed_as_feedback == 1) %>% filter(Target != 'baseline') %>% filter(Muscle=="GM")
succes_rate_GM <- df_subset$Time_percent_success_rate
mean(succes_rate_GM, na.rm=TRUE)
sd(succes_rate_GM, na.rm=TRUE)

df_subset <- df %>% filter(Displayed_as_feedback == 1) %>% filter(Target != 'baseline') %>% filter(Muscle=="GM") %>%
  filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Current_target_imply_size_principle_violation=='FALSE')
succes_rate_GM_compatible_size_principle <- df_subset$Time_percent_success_rate
mean(succes_rate_GM_compatible_size_principle, na.rm=TRUE)
sd(succes_rate_GM_compatible_size_principle, na.rm=TRUE)

df_subset <- df %>% filter(Displayed_as_feedback == 1) %>% filter(Target != 'baseline') %>% filter(Muscle=="GM") %>%
  filter(RT_diff_betwwen_MUs >= min_RT_diff) %>% filter(Current_target_imply_size_principle_violation=='TRUE')
succes_rate_GM_uncompatible_size_principle <- df_subset$Time_percent_success_rate
mean(succes_rate_GM_uncompatible_size_principle, na.rm=TRUE)
sd(succes_rate_GM_uncompatible_size_principle, na.rm=TRUE)

lm_sucess_rate.model <- lmer(Time_percent_success_rate ~ Current_target_imply_size_principle_violation
                             + (Muscle | Subject),
                             data =  df %>% filter(Target != 'baseline') %>% filter(RT_diff_betwwen_MUs >= min_RT_diff))
summary(lm_sucess_rate.model)
lm_sucess_rate.p_val <- anova(lm_sucess_rate.model)
emmeans_obj <- emmeans(lm_sucess_rate.model, ~ Muscle)
pwc_muscles_succes_rate <- pairs(emmeans_obj, adjust = "bonf")
lm_sucess_rate.effsize <- eff_size(emmeans_obj, sigma = sigma(lm_sucess_rate.model),
                                   edf = df.residual(lm_sucess_rate.model))

df_only_GM <- df %>% filter(Muscle=="GM")
df_only_VL <- df %>% filter(Muscle=="VL")

plot_success_rate <- ggplot(data = df %>% filter(Target != 'baseline'),
                            aes(x = Current_target_imply_size_principle_violation,
                                y = Time_percent_success_rate,
                                color = Muscle, fill = Muscle)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               alpha = 0.4, linewidth = 1.5, width = 0.5) +
  
  geom_jitter(data = df_only_GM %>% filter(Target != 'baseline'),
              aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                                              x = -0.15, nudge.from = c("jittered"), seed = 123),
              size = 5, alpha = 0.6) +
  geom_path(data = df_only_GM %>% filter(Target != 'baseline'),
            aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
                                            x = -0.15, nudge.from = c("jittered"), seed = 123),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  geom_jitter(data = df_only_VL %>% filter(Target != 'baseline'),
              aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                                              x = 0.15, nudge.from = c("jittered"), seed = 123),
              size = 5, alpha = 0.6) +
  geom_path(data = df_only_VL %>% filter(Target != 'baseline'),
            aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
                                            x = 0.15, nudge.from = c("jittered"), seed = 123),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  
  scale_color_manual(values=c("dodgerblue","tomato")) +
  scale_fill_manual(values=c("dodgerblue","tomato")) +
  labs(y = "% Success rate", x = "Current target implies size principle violation",
       font = "Roboto",
       title = "Success rate") +
  theme_classic2() + 
  ylim(c(0,100)) +
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
plot_success_rate
ggsave(paste(folder,"Succes_rate_no_RT_diff_threshold.png",sep=""))


#### With RT diff threshold
df_subset <- df %>% filter(RT_diff_betwwen_MUs >= min_RT_diff)
df_only_GM <- df_subset %>% filter(Muscle=="GM")
df_only_VL <- df_subset %>% filter(Muscle=="VL")

lm_sucess_rate.model <- lmer(Time_percent_success_rate ~ Current_target_imply_size_principle_violation
                             + (Muscle | Subject),
                             data =  df_subset %>% filter(Target != 'baseline'))
summary(lm_sucess_rate.model)
lm_sucess_rate.p_val <- anova(lm_sucess_rate.model)
emmeans_obj <- emmeans(lm_sucess_rate.model, ~ Muscle)
pwc_muscles_succes_rate <- pairs(emmeans_obj, adjust = "bonf")
lm_sucess_rate.effsize <- eff_size(emmeans_obj, sigma = sigma(lm_sucess_rate.model),
                                   edf = df.residual(lm_sucess_rate.model))

plot_success_rate_with_RT_diff <- ggplot(data = df_subset %>% filter(Target != 'baseline'),
                            aes(x = Current_target_imply_size_principle_violation,
                                y = Time_percent_success_rate,
                                color = Muscle, fill = Muscle)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               alpha = 0.4, linewidth = 1.5, width = 0.5) +
  
  geom_jitter(data = df_only_GM %>% filter(Target != 'baseline'),
              aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                                              x = -0.15, nudge.from = c("jittered"), seed = 123),
              size = 5, alpha = 0.6) +
  geom_path(data = df_only_GM %>% filter(Target != 'baseline'),
            aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
                                            x = -0.15, nudge.from = c("jittered"), seed = 123),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  geom_jitter(data = df_only_VL %>% filter(Target != 'baseline'),
              aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                  group = MUpair_unique_idx, color = Muscle),
              position = position_jitternudge(width = 0.1,
                                              x = 0.15, nudge.from = c("jittered"), seed = 123),
              size = 5, alpha = 0.6) +
  geom_path(data = df_only_VL %>% filter(Target != 'baseline'),
            aes(x = Current_target_imply_size_principle_violation, y = Time_percent_success_rate,
                group = MUpair_unique_idx, color = Muscle),
            position = position_jitternudge(width = 0.1,
                                            x = 0.15, nudge.from = c("jittered"), seed = 123),
            lineend = "round", linewidth = 1.5, alpha = 0.3) +
  
  
  scale_color_manual(values=c("dodgerblue","tomato")) +
  scale_fill_manual(values=c("dodgerblue","tomato")) +
  labs(y = "% Success rate", x = "Current target implies size principle violation",
       font = "Roboto",
       title = "Success rate, only for (displayed) MU pairs whose RT diff is > X% MVC") +
  theme_classic2() + 
  ylim(c(0,100)) +
  theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
plot_success_rate_with_RT_diff
ggsave(paste(folder,"Succes_rate_with_RT_diff_threshold_3_percent.png",sep=""))


