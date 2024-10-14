# Loading many packages
# install.packages("naniar")
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
library(naniar)

folder <- c("D:/THESE/_DONNEES_MANIP/VL/Number_of_factors_to_select/")
filename <- c("output_table_ALL_SUBJECTS.csv")

df <- read.csv(paste(folder,filename,sep=""),
                      header = TRUE, sep = ",")
#
df$Subject <- factor(df$Subject)
df$Muscle <- factor(df$Muscle)

muscle_levels <- levels(df$Muscle)

#### Rsquared values

temp_df <- df %>% filter(Muscle=="VL")
mean_to_display <- mean(temp_df$Rsquared_associated_below_surrogate)
std_to_display <- sd(temp_df$Rsquared_associated_below_surrogate)
min_to_display <- min(temp_df$Rsquared_associated_below_surrogate)
max_to_display <- max(temp_df$Rsquared_associated_below_surrogate)
print(paste("VL - R² associated with 3 factors = ",round(mean_to_display,2),
            "+-",round(std_to_display,2)," (range: ",
            round(min_to_display,2),"-",round(max_to_display,2),")",
            sep=""))

temp_df <- df %>% filter(Muscle=="GM") %>% filter(Nb_factors_below_surrogate==1)
mean_to_display <- mean(temp_df$Rsquared_associated_below_surrogate)
std_to_display <- sd(temp_df$Rsquared_associated_below_surrogate)
min_to_display <- min(temp_df$Rsquared_associated_below_surrogate)
max_to_display <- max(temp_df$Rsquared_associated_below_surrogate)
print(paste("GM - R² associated with 1 factor = ",round(mean_to_display,2),
            "+-",round(std_to_display,2)," (range: ",
            round(min_to_display,2),"-",round(max_to_display,2),")",
            sep=""))
temp_df <- df %>% filter(Muscle=="GM") %>% filter(Nb_factors_below_surrogate==2)
mean_to_display <- mean(temp_df$Rsquared_associated_below_surrogate)
std_to_display <- sd(temp_df$Rsquared_associated_below_surrogate)
min_to_display <- min(temp_df$Rsquared_associated_below_surrogate)
max_to_display <- max(temp_df$Rsquared_associated_below_surrogate)
print(paste("GM - R² associated with 2 factors = ",round(mean_to_display,2),
            "+-",round(std_to_display,2)," (range: ",
            round(min_to_display,2),"-",round(max_to_display,2),")",
            sep=""))

### Number of factors

temp_df <- df %>% filter(Muscle=="VL")
values_to_use <- temp_df$Nb_factors_5percent_threshold
values_to_use <- append(values_to_use, temp_df$Nb_factors_linear_fit)
mean_to_display <- mean(values_to_use)
std_to_display <- sd(values_to_use)
min_to_display <- min(values_to_use)
max_to_display <- max(values_to_use)
print(paste("VL - Number of factors found by the two other methods = ",round(mean_to_display,2),
            "+-",round(std_to_display,2)," (range: ",
            round(min_to_display,2),"-",round(max_to_display,2),")",
            sep=""))

temp_df <- df %>% filter(Muscle=="GM")
values_to_use <- temp_df$Nb_factors_5percent_threshold
values_to_use <- append(values_to_use, temp_df$Nb_factors_linear_fit)
mean_to_display <- mean(values_to_use)
std_to_display <- sd(values_to_use)
min_to_display <- min(values_to_use)
max_to_display <- max(values_to_use)
print(paste("GM - Number of factors found by the two other methods = ",round(mean_to_display,2),
            "+-",round(std_to_display,2)," (range: ",
            round(min_to_display,2),"-",round(max_to_display,2),")",
            sep=""))

######################## PLOT R² FACTOR ANALYSIS (true and surrogate data)
# # # # # # # 
folder <- c("D:/THESE/_DONNEES_MANIP/VL/Number_of_factors_to_select/")
filename_true_data <- c("Rsquared_for_different_factor_nb_ALL_SUBJECTS.csv")
filename_surrogate_data <- c("Rsquared_for_different_factor_nb_surrogate_ALL_SUBJECTS.csv")


df_experiment <- read.csv(paste(folder,filename_true_data,sep=""),
               header = TRUE, sep = ",")
df_experiment$data <- "experimental"
df_surrogate <- read.csv(paste(folder,filename_surrogate_data,sep=""),
                          header = TRUE, sep = ",")
df_surrogate$data <- "surrogate"

#
df_all_factors <- rbind(df_experiment, df_surrogate)
df_all_factors$data <- factor(df_all_factors$data)
df_all_factors$Subject <- factor(df_all_factors$Subject)
df_all_factors$Muscle <- factor(df_all_factors$Muscle)

muscle_levels <- levels(df_all_factors$Muscle)

### Reorganize table so that each nb of factor is in one variable
df_new <- data.frame(Subject = factor(),
                     Muscle = factor(),
                     Nb_factors = numeric(),
                     Rsquared = numeric(),
                     data = factor()
)
for (subjecti in 1:nrow(df_all_factors)) {
  for (factori in 3:ncol(df_all_factors)) {
    new_row <- data.frame(
      Subject = df_all_factors$Subject[subjecti],
      Muscle = df_all_factors$Muscle[subjecti],
      Nb_factors = factori-3,
      Rsquared = df_all_factors[subjecti,factori],
      data = df_all_factors$data[subjecti]
    )
    df_new <- rbind(df_new,new_row)
  }
}
df_new$Rsquared <- as.numeric(df_new$Rsquared)

#### PLOT R² according to nb of factors
df_new <- df_new %>% filter(Muscle=="VL") %>% #"GM"
  filter(Nb_factors < 6)
# Reorder the levels of Nb_factors based on the mean Rsquared values
mean_Rsquared <-  df_new %>% filter(Nb_factors >= 1) %>% filter(Nb_factors < 3) %>%
  filter(data == "experimental") %>% group_by(Subject) %>%
  summarize(mean_Rsquared = mean(Rsquared, na.rm = TRUE))
mean_Rsquared$Subject <- droplevels(mean_Rsquared$Subject)
df_new$Subject <- droplevels(df_new$Subject)
levels(df_new$Subject)
levels(mean_Rsquared$Subject)
mean_Rsquared$Subject <- factor(mean_Rsquared$Subject, 
                         levels = mean_Rsquared$Subject[order(mean_Rsquared$mean_Rsquared)])
df_new$Subject <- factor(df_new$Subject, 
                            levels = levels(mean_Rsquared$Subject))

data_to_plot <- df_new %>% filter(data == "experimental")
data_to_plot_surrogate <- df_new %>% filter(data == "surrogate")
data_nb_factors_to_plot <- df %>% filter(Muscle == "VL") #"GM"

# Create offset for each Subject
dodge_per_subject <- c((1:length(levels(data_to_plot$Subject)))-
                         (length(levels(data_to_plot$Subject))+1)/2)*
                          0
for (rowi in 1:nrow(data_to_plot)) {
  subject_idex <- which(levels(data_to_plot$Subject) == data_to_plot$Subject[rowi])
  data_to_plot$Nb_factors[rowi] = data_to_plot$Nb_factors[rowi] + dodge_per_subject[subject_idex]
  data_to_plot_surrogate$Nb_factors[rowi] = data_to_plot_surrogate$Nb_factors[rowi] + dodge_per_subject[subject_idex]
}
for (rowi in 1:nrow(data_nb_factors_to_plot)) {
  subject_idex <- which(levels(data_to_plot$Subject) == data_nb_factors_to_plot$Subject[rowi])
  data_nb_factors_to_plot$Nb_factors_below_surrogate[rowi] = data_nb_factors_to_plot$Nb_factors_below_surrogate[rowi] + dodge_per_subject[subject_idex]
}
# update offset for data_nb_of_factors_to_plot

# Define a set of custom colors
custom_colors <- c("tomato","firebrick4") # c("dodgerblue", "royalblue4")
custom_fills <- c("tomato","firebrick4") # c("dodgerblue", "royalblue4")
# custom_colors <- c("tomato", "tomato4")
# Create the color ramp palette function
custom_palette_function <- colorRampPalette(custom_colors)
custom_colormap <- custom_palette_function(length(levels(data_to_plot$Subject))) # don't know why I need +1 but otherwise I get an error
custom_palette_function <- colorRampPalette(custom_fills)
custom_fillmap <- custom_palette_function(length(levels(data_to_plot$Subject))) # don't know why I need +1 but otherwise I get an error

plot_Rsquared <- ggplot(data = data_to_plot,
                        aes(x = Nb_factors, y = Rsquared,  color = Subject,
                            fill = Subject, group = Subject)) +
  geom_line(data = data_to_plot_surrogate,
            aes(x = Nb_factors, y = Rsquared,  color = Subject,
                group = Subject), linetype = "longdash",
            linewidth = 1.2, alpha = 1) + # surrogate data lines
  geom_line(linewidth = 1.2, alpha = 1) +
  geom_point(size = 5, alpha = 1, shape = 21, stroke = 2.5) +
  #geom_point(size = 6, alpha = 1, shape = 21, stroke = 2.5, fill=NA) +
  geom_point(data = data_nb_factors_to_plot,
             aes(x = Nb_factors_below_surrogate, y = Rsquared_associated_below_surrogate, 
                 color = Subject),
             size = 6, alpha = 1, shape = 21,
             stroke = 3, fill = "white") +
  scale_color_manual(values=custom_colormap) +
  scale_fill_manual(values=custom_fillmap) +
  labs(x = "Number of latent factors", y = "R²", font = "Poppins") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 1)) +  # Hard cutoff for axes
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1.5, color = "darkgrey"),
    axis.title.x = element_text(size = 20, color = "grey35", face = "bold"), # Axis titles
    axis.title.y = element_text(size = 20, color = "grey35", face = "bold",
                                angle = 0, vjust = 0.5), # Axis titles
    axis.text.x = element_text(size = 20, color = "grey40", face = "bold"), # Axis text
    axis.text.y = element_text(size = 20, color = "grey40", face = "bold"), # Axis text
    axis.ticks = element_line(color = "darkgrey", linewidth = 1.5), # Customize tick color and size
    axis.ticks.length = unit(0.25,"cm"),
    legend.position = "none"
  )

plot_Rsquared
ggsave(filename = "plot_Rsquared_GM.png", plot = plot_Rsquared, 
       width = 2300, height = 2300, units = "px")
# ggsave(filename = "plot_Rsquared_GM.png", plot = plot_Rsquared_GM, 
#        width = 3000, height = 2863, units = "px")


means_data_to_plot <- aggregate(
  Rsquared ~ Nb_factors + Muscle, data = data_to_plot, FUN = function(x) c(mean = mean(x), std = sd(x)))

# plot_Rsquared <- ggplot(data = data_to_plot,
#                         aes(x = Nb_factors, y = Rsquared, color = Muscle)) +
#   geom_line(data = means_data_to_plot,
#             aes(x = Nb_factors, y = Rsquared[,"mean"], color = Muscle),
#             linewidth = 2) +
#   geom_errorbar(data = means_data_to_plot,
#                 aes(x = Nb_factors, y = Rsquared[,"mean"],
#                     ymin = Rsquared[,"mean"] - Rsquared[,"std"],
#                     ymax = Rsquared[,"mean"] + Rsquared[,"std"],
#                     color = Muscle),
#                 linewidth = 1.5, width = 0.2, alpha = 0.3) +
#   geom_jitter(size = 4, alpha = 0.6, width = 0.05) +
#   scale_color_manual(values=c("dodgerblue","tomato")) +
#   labs(x = "Number of factors", y = "R²", font = "Poppins") +
#   theme_classic2() + 
#   theme(axis.line = element_line(linewidth = 1, color = "darkgrey"))
# plot_Rsquared
# ggsave(paste(folder_save,"Rsquared_per_nb_factors_RAW.png",sep=""))
