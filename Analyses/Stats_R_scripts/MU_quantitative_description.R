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

df <- read.csv("D:/THESE/_ARTICLES_WRITING/General_descriptive_report.csv",
         header = TRUE, sep = ",")

df$Subject <- as.factor(df$Subject)
df$Muscle <- as.factor(df$Muscle)
df$Condition <- as.factor(df$Condition)

### TOTAL individual MUs decomposed across conditions and participants
sum(df$total_independent_MUs_decoded_across_conditions, na.rm = TRUE)
### TOTAL individual MUs tracked at least once, across conditions and participants
sum(df$total_MUs_matched_at_least_once, na.rm = TRUE)

print("NUMBER OF MUs USED FOR FA")
for (musclei in 1:2)
{
  if (musclei==1)
  {
    current_mucle_considered <- "VL"
  } else {
    current_mucle_considered <- "GM"
  }
  df_to_use <- df %>% filter(Muscle==current_mucle_considered)
  df_to_use_temp <- df_to_use %>% filter(Condition==levels(df_to_use$Condition)[1]) # only plateau
  
  print("NUMBER ------ ")
  print(paste("Mean # of MUs used for FA for ",current_mucle_considered,
              " = ",round(mean(df_to_use_temp$MUs_continuous_on_ref),2),
              sep=""))
  print(paste("StD # of MUs used for FA for ",current_mucle_considered,
              " = ",round(sd(df_to_use_temp$MUs_continuous_on_ref),2),
              sep=""))
  print(paste("Min # of MUs used for FA for ",current_mucle_considered,
              " = ",round(min(df_to_use_temp$MUs_continuous_on_ref),2),
              sep=""))
  print(paste("Max # of MUs used for FA for ",current_mucle_considered,
              " = ",round(max(df_to_use_temp$MUs_continuous_on_ref),2),
              sep=""))
  print("")
  print("RATIO ------ ")
  print(paste("Mean % of MUs used for FA for ",current_mucle_considered,
              " = ",round(mean(df_to_use_temp$MUs_continuous_on_ref/df_to_use_temp$MUs_decoded)*100,1),"%",
              sep=""))
  print(paste("StD % of MUs used for FA for ",current_mucle_considered,
              " = ",round(sd(df_to_use_temp$MUs_continuous_on_ref/df_to_use_temp$MUs_decoded)*100,1),"%",
              sep=""))
  print(paste("Min % of MUs used for FA for ",current_mucle_considered,
              " = ",round(min(df_to_use_temp$MUs_continuous_on_ref/df_to_use_temp$MUs_decoded)*100,1),"%",
              sep=""))
  print(paste("Max % of MUs used for FA for ",current_mucle_considered,
              " = ",round(max(df_to_use_temp$MUs_continuous_on_ref/df_to_use_temp$MUs_decoded)*100,1),"%",
              sep=""))
  print("")
}

current_mucle_considered <- "VL"
df_to_use <- df %>% filter(Muscle==current_mucle_considered)

for (conditioni in 1:(length(levels(df_to_use$Condition))-1) )
{
  df_to_use_temp <- df_to_use %>% filter(Condition==levels(df_to_use$Condition)[conditioni])
  print(paste("Mean for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",round(mean(df_to_use_temp$MUs_decoded),2),
              sep=""))
  print(paste("StD for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",round(sd(df_to_use_temp$MUs_decoded),2),
              sep=""))
  print(paste("Minimum for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",min(df_to_use_temp$MUs_decoded),
              sep=""))
  print(paste("Maximum for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",max(df_to_use_temp$MUs_decoded),
              sep=""))
  print("")
}

for (conditioni in 1:(length(levels(df_to_use$Condition))-1) )
{
  df_to_use_temp <- df_to_use %>% filter(Condition==levels(df_to_use$Condition)[conditioni])
  print(paste("Mean for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",round(mean(df_to_use_temp$MUs_matched_with_ref),2),
              sep=""))
  print(paste("StD for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",round(sd(df_to_use_temp$MUs_matched_with_ref),2),
              sep=""))
  print(paste("Minimum for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",min(df_to_use_temp$MUs_matched_with_ref),
              sep=""))
  print(paste("Maximum for ",current_mucle_considered,
              " condition ",levels(df_to_use$Condition)[conditioni],
              " = ",max(df_to_use_temp$MUs_matched_with_ref),
              sep=""))
  print("")
}


# total decoded
for (muscli in 1:length(levels(df$Muscle)))
{
  df_to_use_temp <- df %>% filter(Muscle==levels(Muscle)[muscli]) %>% filter(Condition==levels(df_to_use$Condition)[length(levels(df_to_use$Condition))])
  current_mucle_considered = levels(df$Muscle)[muscli]
  print(paste("Mean for ",current_mucle_considered,
              " = ",round(mean(df_to_use_temp$total_independent_MUs_decoded_across_conditions),2),
              sep=""))
  print(paste("StD for ",current_mucle_considered,
              " = ",round(sd(df_to_use_temp$total_independent_MUs_decoded_across_conditions),2),
              sep=""))
  print(paste("Minimum for ",current_mucle_considered,
              " = ",min(df_to_use_temp$total_independent_MUs_decoded_across_conditions),
              sep=""))
  print(paste("Maximum for ",current_mucle_considered,
              " = ",max(df_to_use_temp$total_independent_MUs_decoded_across_conditions),
              sep=""))
  print("")
  df_to_use_temp <- df %>% filter(Muscle==levels(Muscle)[muscli]) %>% filter(Condition=="total")
  print(paste("Total number of MUs matched at least once for ",current_mucle_considered,
              " across subjects and conditions = ",sum(df_to_use_temp$total_MUs_matched_at_least_once),
              sep=""))
  print(paste("Mean number of MUs matched at least once for ",current_mucle_considered,
              " across subjects and conditions = ",round(
                mean(df_to_use_temp$total_MUs_matched_at_least_once),0),
              " +- ", round(sd(df_to_use_temp$total_MUs_matched_at_least_once),0),
              sep=""))
  print("")
  print("")
  print("")
}

# MUs matched in all files
for (muscli in 1:length(levels(df$Muscle)))
{
  df_to_use_temp <- df %>% filter(Muscle==levels(Muscle)[muscli]) %>% filter(Condition==levels(df_to_use$Condition)[length(levels(df_to_use$Condition))])
  current_mucle_considered = levels(df$Muscle)[muscli]
  print(paste("Mean for ",current_mucle_considered,
              " = ",round(mean(df_to_use_temp$total_MUs_matched_in_all_conditions),2),
              sep=""))
  print(paste("StD for ",current_mucle_considered,
              " = ",round(sd(df_to_use_temp$total_MUs_matched_in_all_conditions),2),
              sep=""))
  print(paste("Minimum for ",current_mucle_considered,
              " = ",min(df_to_use_temp$total_MUs_matched_in_all_conditions),
              sep=""))
  print(paste("Maximum for ",current_mucle_considered,
              " = ",max(df_to_use_temp$total_MUs_matched_in_all_conditions),
              sep=""))
  print("")
}

### % matched relative to plateau condition
# MUs matched in all files
print("% MUs matched relative to the plateau condition (on which FA was performed)")
for (conditioni in 2:(length(levels(df$Condition))-1) ) # skip plateau and total
{
  print(paste("Condition ",levels(df$Condition)[conditioni]," =============",sep=""))
  print("")
  for (muscli in 1:length(levels(df$Muscle))) 
  {
    # only VL and GM
    if ((muscli == 2) || (muscli == 4))
    {
      next
    }
    df_to_use_temp <- df %>% filter(Muscle==levels(Muscle)[muscli]) %>% filter(Condition==levels(df_to_use$Condition)[conditioni])
    current_mucle_considered = levels(df$Muscle)[muscli]
    ref_to_use_temp <- df %>% filter(Muscle==levels(Muscle)[muscli]) %>% filter(Condition==levels(df_to_use$Condition)[1])
    # print(paste("Mean for ",current_mucle_considered,
    #             " = ",round(mean(df_to_use_temp$MUs_matched_with_ref),2),
    #             sep=""))
    # print(paste("StD for ",current_mucle_considered,
    #             " = ",round(sd(df_to_use_temp$MUs_matched_with_ref),2),
    #             sep=""))
    # print(paste("Minimum for ",current_mucle_considered,
    #             " = ",min(df_to_use_temp$MUs_matched_with_ref),
    #             sep=""))
    # print(paste("Maximum for ",current_mucle_considered,
    #             " = ",max(df_to_use_temp$MUs_matched_with_ref),
    #             sep=""))
    # print("")
    print(paste("Mean for ",current_mucle_considered,
                " = ",round(
                  mean(df_to_use_temp$MUs_matched_with_ref /
                    ref_to_use_temp$MUs_decoded * 100)
                    ,2),"%", sep=""))
    print(paste("StD for ",current_mucle_considered,
                " = ",round(
                  sd(df_to_use_temp$MUs_matched_with_ref /
                         ref_to_use_temp$MUs_decoded * 100)
                  ,2),"%", sep=""))
    print(paste("Minimum for ",current_mucle_considered,
                " = ",round(
                  min(df_to_use_temp$MUs_matched_with_ref /
                         ref_to_use_temp$MUs_decoded * 100)
                  ,2),"%", sep=""))
    print(paste("Maximum for ",current_mucle_considered,
                " = ",round(
                  max(df_to_use_temp$MUs_matched_with_ref /
                         ref_to_use_temp$MUs_decoded * 100)
                  ,2),"%", sep=""))
    print("")
  }
  print("")
  print("")
}
