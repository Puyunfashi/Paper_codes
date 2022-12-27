# To conduct two-way ANOVA on the scale scores of affects after each session

# Clear environment
rm(list = ls()) 

# Set system language
Sys.setenv(LANG = "en")

# LOAD PACKAGES ############################################
# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Then load the package by using either of the following:
library(pacman)  # No message.

# Or, by using "pacman::p_load" you can use the p_load
# function from pacman without actually loading pacman.
# These are packages I load every time.
pacman::p_load(pacman, GGally, ggplot2, ggthemes, ggsignif, Hmisc,
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, ggwordcloud, svglite, rstatix, dplyr, tidyverse) 

work_dir <- "/Users/ChenXiao/Documents/My_Documents/rumination_coorp/Suzhou_analysis/behavioral_results"
rio_csv <- import(paste(work_dir, "/behaviral_results.csv", sep = "")) 

rio_csv_NA <- rio_csv %>%                               # Replacing values
  mutate(Rest_emo_before = replace(Rest_emo_before, Rest_emo_before == 999, NA),
         Rest_emo_after = replace(Rest_emo_after, Rest_emo_after == 999, NA),
         Sad_emo = replace(Sad_emo, Sad_emo == 999, NA),
         Rum_emo = replace(Rum_emo, Rum_emo == 999, NA),
         Dis_emo = replace(Dis_emo, Dis_emo == 999, NA))

rio_csv_noNA <-  rio_csv_NA %>%
  filter(!is.na(Rest_emo_before) & !is.na(Rest_emo_after) &
           !is.na(Sad_emo) & !is.na(Rum_emo) & !is.na(Dis_emo))

rio_csv_noNA_long <- rio_csv_noNA %>% 
  pivot_longer(
    cols=c('Rest_emo_before', 'Rest_emo_after', 'Sad_emo', 'Rum_emo', 'Dis_emo'),
    names_to = "Condition",
    values_to = "Score"
  )

rio_csv_noNA_long <- rio_csv_noNA_long %>% 
convert_as_factor(Dx, Condition)
rio_csv_noNA_long <- rio_csv_noNA_long %>%
reorder_levels("Condition", order = c("Rest_emo_before", "Rest_emo_after",
                                      "Sad_emo", "Rum_emo", "Dis_emo"))
                                       
################ do statistics ##################
res.aov <- anova_test(
  data = rio_csv_noNA_long, dv = Score, wid = SubList,
  within =  Condition,
  between = Dx
)
get_anova_table(res.aov)

# Effect of treatment at each time point
one.way <- rio_csv_noNA_long %>%
  group_by(Dx) %>%
  anova_test(dv = Score, wid = SubList, within = Condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

# independent comparisons between Dx groups
pwc <- rio_csv_noNA_long %>%
  group_by(Condition) %>%
  t_test(
    Score ~ Dx, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

# Effect of time at each level of treatment
one.way2 <- rio_csv_noNA_long %>%
  group_by(Dx) %>%
  anova_test(dv = Score, wid = SubList, within = Condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2
# Pairwise comparisons between time points
pwc2 <- rio_csv_noNA_long %>%
  group_by(Dx) %>%
  pairwise_t_test(
    Score ~ Condition, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2

# plot figure
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df4plot <- rio_csv_noNA_long
df4plot$Condition <- recode_factor(df4plot$Condition, 
                                      Rest_emo_before = "Before_Resting", 
                                      Rest_emo_after = "After_Resting",
                                      Sad_emo = "Sad", Rum_emo = "Rum", Dis_emo = "Dis")
colnames(df4plot)[2] <- "Group"
df4plot$Group <- recode_factor(df4plot$Group, "1" = "MDD", "2" = "HC")

df3 <- data_summary(df4plot, varname="Score", 
                    groupnames=c("Group", "Condition"))
head(df3)

# dot and line
ggplot(data = df3, aes(Condition, Score, color = Group, fill = Group)) +
  theme_classic(base_size = 20) + 
  geom_errorbar(aes(ymin=Score-sd, ymax=Score+sd), width= .5, size = 1,  
                position=position_dodge(0.05)) +
  geom_line(aes(linetype=Group, group = Group), size = 1) + 
  geom_point(aes(shape=Group), size = 4) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 12)) +
  geom_signif(comparisons = list(c("Before_Resting", "Sad"),
                                 c("Before_Resting", "Rum")), 
              annotations=c("***", "*"),
              map_signif_level = FALSE, textsize = 6, color = "black")

# violin plot
ggplot(data = df4plot, aes(Condition, Score, color = Group, fill = Group)) +
  theme_classic(base_size = 32) + xlab("") +
  geom_violin(aes(Condition, Score, color = Group, fill = Group)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 32)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
               geom="pointrange", color="black",
               shape = 18, size = 2,
               position = position_dodge(width = 0.9))

ggsave(paste(work_dir, "/Emotion_ratings.png", sep = ""), 
       width = 30, height = 20, units = "cm", dpi = 300)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear packages
p_unload(all)  # Easier: clears all add-ons
detach("package:datasets", unload = TRUE)  # For base
