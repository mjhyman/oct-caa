# load necessary libraries
library(nlme)
library(dplyr)
library(readxl)
library(ggplot2)
library(rstanarm)
library(tidybayes)
library(tidyverse)
library(lme4)

### SECTION 3: EDA ###

# # original retardance data from October
# 
# retardance_data_eda <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "retardance")
# 
# # rename OpticalProperty to be retardance
# retardance_data_eda <- retardance_data_eda %>%                                                    
#   rename(retardance = OpticalProperty)
# 
# # Figure 1: side by side boxplots of retardance by group
# ggplot(retardance_data_eda, aes(x = Groups, y = log(retardance), fill=Groups)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "Log-Transformed Retardance by Group", x = "Group", y = "log(Retardance)") 
# 
# summary(retardance_data_eda$retardance)
# 
# # original scattering coefficient data from October
# scattering_data_eda <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "scattering")
# 
# # rename OpticalProperty to be scattering coefficient
# scattering_data_eda <- scattering_data_eda %>%                                                    
#   rename(scattering = OpticalProperty) 
# 
# # Figure 2: side by side boxplots of scattering coefficient by group
# ggplot(scattering_data_eda, aes(x = Groups, y = log(scattering), fill=Groups)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "Log-Transformed Scattering Coefficient by Group", x = "Group", y = "log(Scattering Coefficient)") 
# 
# summary(scattering_data_eda$scattering)
# 
# ### The rest of the report uses updated data from November ###
# # EDA Updated Data
# # load in necessary data sets
# 
# # retardance data set
# retardance_data <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_03Nov2025.xlsx", sheet = "retardance")
# 
# # Figure 3
# ggplot(retardance_data, aes(x = Region, fill = Groups)) +
#   geom_bar(position = position_dodge()) +
#   geom_text(
#     stat = "count",
#     aes(label = after_stat(count)),
#     position = position_dodge(width = 0.9),
#     vjust = -0.3,
#     size = 3
#   ) +
#   labs(
#     title = "Number of Observations for Retardance",
#     x = "Region",
#     y = "Count of Observations"
#   ) +
#   theme_minimal()
# 
# # Read in the scattering dataset
# scattering_data <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_03Nov2025.xlsx", sheet = "scattering")
# 
# # Figure 4
# ggplot(scattering_data, aes(x = Region, fill = Groups)) +
#   geom_bar(position = position_dodge()) +
#   geom_text(
#     stat = "count",
#     aes(label = after_stat(count)),
#     position = position_dodge(width = 0.9),
#     vjust = -0.3,
#     size = 3
#   ) +
#   labs(
#     title = "Number of Observations for Scattering Coefficient",
#     x = "Region",
#     y = "Count of Observations"
#   ) +
#   theme_minimal()
# 
# # Figure 5
# count_retardance_data <- retardance_data %>%
#   count(distance, Region)
# 
# 
# ggplot(count_retardance_data,
#        aes(x = distance, y = n,
#            color = Region, group = Region)) +
#   geom_line() +
#   labs(
#     title = "Sample Size across Distance by Region",
#     x = "Distance",
#     y = "Count (log scale)",
#     color = "Region"
#   ) +
#   scale_y_log10() +
#   theme_minimal()

# Figure 6
#retardance_long <- retardance_data %>%
#  mutate(property = "Retardance") %>%
#  rename(value = retardance)   # rename the measurement to a common column

#scattering_long <- scattering_data %>%
#  mutate(property = "Scattering") %>%
#  rename(value = scattering)

# Combine datasets
# combined_data <- bind_rows(scattering_long, retardance_long)
# 
# ggplot(combined_data, aes(x = Groups, y = value, fill = Groups)) +
#   geom_boxplot() +
#   facet_grid(property ~ Region, scales = "free_y") +
#   labs(
#     title = "Side-by-side Boxplots by Group, Region, and Optical Property",
#     x = "Group",
#     y = "log(Optical Property)",
#     fill = "Group"
#   ) +
#   theme_minimal() 
# 
# # Figure 7
# scattering_data_mean <- scattering_data %>%
#   group_by(Groups, distance) %>%
#   summarise(mean_opt = mean(scattering), .groups = "drop")
# 
# ggplot(scattering_data_mean, aes(x = distance, y = mean_opt, color = Groups)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Mean optical property vs distance",
#        x = "Distance from vessel/EPVS (￿m)",
#        y = "Mean optical property",
#        color = "Group")
# 
# # Figure 8
# retardance_data_mean <- retardance_data %>%
#   group_by(Groups, distance) %>%
#   summarise(mean_opt = mean(retardance), .groups = "drop")
# 
# ggplot(retardance_data_mean, aes(x = distance, y = mean_opt, color = Groups)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Mean optical property vs distance",
#        x = "Distance from vessel/EPVS (￿m)",
#        y = "Mean optical property",
#        color = "Group")
# 
# 
# # Figure 9
# caa_mean_region <- scattering_data %>%
#   group_by(Region, Groups, distance) %>%
#   summarise(mean_opt = mean(scattering), .groups = "drop")
# 
# ggplot(caa_mean_region,
#        aes(x = distance, y = mean_opt, color = Groups)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~ Region) +
#   labs(title = "Mean scattering coefficient vs distance by region",
#        x = "Distance from vessel/EPVS (￿m)",
#        y = "Mean optical property",
#        color = "Group")
# 
# # Figure 10
# caa_mean_region <- retardance_data %>%
#   group_by(Region, Groups, distance) %>%
#   summarise(mean_opt = mean(retardance), .groups = "drop")
# 
# ggplot(caa_mean_region,
#        aes(x = distance, y = mean_opt, color = Groups)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~ Region) +
#   labs(title = "Mean optical property vs distance by region",
#        x = "Distance from vessel/EPVS (￿m)",
#        y = "Mean optical property",
#        color = "Group")








### SECTION 4: STATISTICAL ANALYSIS ###
### Part 1 ###
### load in necessary data sets ###

# retardance data set
retardance_data <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_03Nov2025.xlsx", sheet = "retardance")

# rename OpticalProperty to be retardance
retardance_data <- retardance_data %>%                                                    
  rename(retardance = OpticalProperty)

# dataset with retardance means across Region, Groups, and distance
retardance_means <- retardance_data %>%
  group_by(Region,subID, Groups, distance) %>%
  summarise(
    mean_value = mean(retardance, na.rm = TRUE),
    .groups = "drop"
  )

# Create a factor version of distance
retardance_means <- retardance_means %>%
  mutate(distance_f = factor(distance))
head(retardance_means)


# Read in the scattering dataset
scattering_data <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_03Nov2025.xlsx", sheet = "scattering")

# rename OpticalProperty to be scattering 
scattering_data <- scattering_data %>%                                                    
  rename(scattering = OpticalProperty)



### Approach 1: Frequentist ###

### Stage 1: Linear mixed effect model assuming no spatial correlation
# scattering coefficient frontal region
assump_frontal_data_scattering <- scattering_data %>% 
  filter(Region == "front") %>%
  mutate(distance_f = factor(distance))

mod_front_scat <- lmer(
  scattering ~ Groups * distance_f + (1 | subID),
  data = assump_frontal_data_scattering
)

# table 1 of appendix 
summary(mod_front_scat)

# scattering coefficient occipital region
assump_occip_data_scattering <- scattering_data %>% 
  filter(Region == "occip") %>%
  mutate(distance_f = factor(distance))

mod_occip_scat <- lmer(
  scattering ~ Groups * distance_f + (1 | subID),
  data = assump_occip_data_scattering
)

# table 2 of appendix 
summary(mod_occip_scat)

# scattering coefficient combined region
mod_comb_scat <- lmer(
  scattering ~ Groups * factor(distance) + Region + (1 | subID) + (1 | subID:Region),
  data = scattering_data
)

# table 3 of appendix 
summary(mod_comb_scat)

# retardance frontal region
assump_frontal_data_retardance <- retardance_data %>% 
  filter(Region == "front") %>%
  mutate(distance_f = factor(distance))

mod_front_ret <- lmer(
  retardance ~ Groups * distance_f + (1 | subID),
  data = assump_frontal_data_retardance
)

# table 4 of appendix 
summary(mod_front_ret)

# retardance occipital region
assump_occip_data_retardance <- retardance_data %>% 
  filter(Region == "occip") %>%
  mutate(distance_f = factor(distance))

mod_occip_ret <- lmer(
  retardance ~ Groups * distance_f + (1 | subID),
  data = assump_occip_data_retardance
)

# table 5 of appendix 
summary(mod_occip_ret)

# retardance combined region
mod_comb_ret <- lmer(
  retardance ~ Groups * factor(distance) + Region + (1 | subID) + (1 | subID:Region),
  data = retardance_data
)

# table 6 of appendix 
summary(mod_comb_ret)

# Figure 11
par(mfrow = c(2, 3))
qqnorm(residuals(mod_front_scat),
       main = "QQ: Scattering (Frontal)")
qqline(residuals(mod_front_scat))
qqnorm(residuals(mod_occip_scat),
       main = "QQ: Scattering (Occipital)")
qqline(residuals(mod_occip_scat))
qqnorm(residuals(mod_comb_scat),
       main = "QQ: Scattering (Combined)")
qqline(residuals(mod_front_scat))
qqnorm(residuals(mod_front_ret),
       main = "QQ: Retardance (Frontal)")
qqline(residuals(mod_front_ret))
qqnorm(residuals(mod_occip_ret),
       main = "QQ: Retardance (Occipital)")
qqline(residuals(mod_occip_ret))
qqnorm(residuals(mod_comb_ret),
       main = "QQ: Scattering (Combined)")
qqline(residuals(mod_comb_ret))

### Stage 2: Linear mixed effect model assuming spatial correlation
# This was only done for combined regions

# dataset with retardance means across Region, Groups, and distance
retardance_means <- retardance_data %>%
  group_by(Region,subID, Groups, distance) %>%
  summarise(
    mean_value = mean(retardance, na.rm = TRUE),
    .groups = "drop"
  )

# dataset with scattering coef means across Region, Groups, and distance
scattering_means <- scattering_data %>%
  group_by(Region, subID, Groups, distance) %>%
  summarise(
    mean_value = mean(scattering, na.rm = TRUE),
    .groups = "drop"
  )

# Create a factor version of distance
retardance_means <- retardance_means %>%
  mutate(distance_f = factor(distance))
head(retardance_means)

# Create a factor version of distance
scattering_means <- scattering_means %>%
  mutate(distance_f = factor(distance))
head(scattering_means)

## Combined Retardance ##

# Create a new variable subjectID that corresponds to the 5 subjects
# rather than the regions
# this is my best guess but needs to be checked
m_retardance_data <- retardance_data %>%
  mutate(subjectID = case_when(
    subID %in% c(1, 2) ~ 1,
    subID == 3        ~ 2,
    subID %in% c(4, 5) ~ 3,
    subID %in% c(6, 7) ~ 4,
    subID %in% c(8, 9) ~ 5
  ))

table(m_retardance_data$subjectID, m_retardance_data$Region)

# average retardance across Groups, Region, subjectID,and distance
means_mret_data <- m_retardance_data %>%
  group_by(Groups, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(retardance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

# fit classical linear mixed-effects model with random intercept
mod_comb_ret <- lmer(
  mean_value ~ Groups * Region * distance_f + (1 | subjectID), # random intercept for subjectID 
  data = means_mret_data
)

# Table 12: Check Normality
qqnorm(residuals(mod_comb_ret),
       main = "QQ plot: Combined Retardance (residuals)")
qqline(residuals(mod_comb_ret))


### Approach 2: Bayesian ###
### Models for Combined regions ###

## Combined Retardance ##

# Create a new variable subjectID that corresponds to the 5 subjects
# rather than the regions
# this is my best guess but needs to be checked
m_retardance_data <- retardance_data %>%
  mutate(subjectID = case_when(
    subID %in% c(1, 2) ~ 1,
    subID == 3        ~ 2,
    subID %in% c(4, 5) ~ 3,
    subID %in% c(6, 7) ~ 4,
    subID %in% c(8, 9) ~ 5
  ))

table(m_retardance_data$subjectID, m_retardance_data$Region)

# average retardance across Groups, Region, subjectID,and distance
means_mret_data <- m_retardance_data %>%
  group_by(Groups, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(retardance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

# fit a Bayesian linear mixed-effects model
bmod_comb_ret <- stan_glmer(
  mean_value ~ Groups * Region * distance_f + (1 | subjectID),
  data = means_mret_data,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  cores = 4
)

# Figure 13
pp_check(bmod_comb_ret)

pp_check(bmod_comb_ret, "intervals") +
  facet_wrap(means_mret_data$Groups ~ means_mret_data$Region)

pp_check(bmod_comb_ret, "intervals") +
  facet_wrap(~ means_mret_data$Region)

pp_check(bmod_comb_ret, "intervals") +
  facet_wrap(~ means_mret_data$Groups)

## Combined Scattering ##

# Create a new variable subjectID that corresponds to the 5 subjects
# rather than the regions
# this is my best guess but needs to be checked
m_scattering_data <- scattering_data %>%
  mutate(subjectID = case_when(
    subID %in% c(1, 2) ~ 1,
    subID == 3        ~ 2,
    subID %in% c(4, 5) ~ 3,
    subID %in% c(6, 7) ~ 4,
    subID %in% c(8, 9) ~ 5,
    TRUE ~ NA_real_
  )) 

table(m_scattering_data$subjectID, m_scattering_data$Region)

# average scattering coefficient across Groups, Region, subjectID, and distance
means_mscat_data <- m_scattering_data %>%
  group_by(Groups, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(scattering, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

# Next fit a Bayesian linear mixed-effects model
bmod_comb_scat <- stan_glmer(
  mean_value ~ Groups * Region * distance_f + (1 | subjectID),
  data = means_mscat_data,
  family = gaussian,
  chains = 4,
  iter = 4000,
  cores = 4
)

# Figure 14
pp_check(bmod_comb_scat)

pp_check(bmod_comb_scat, "intervals") +
  facet_wrap(means_mscat_data$Groups ~ means_mscat_data$Region)

pp_check(bmod_comb_scat, "intervals") +
  facet_wrap(~ means_mscat_data$Region)

pp_check(bmod_comb_scat, "intervals") +
  facet_wrap(~ means_mscat_data$Groups)

# Figure 15: Bias per Subject Across Distance Faceted by Group and Region (Retardance)
bias_df_ret <- means_mret_data %>%
  mutate(
    pred_mean = fitted(bmod_comb_ret),       
    bias      = mean_value - pred_mean
  )

# Plot bias vs distance
ggplot(bias_df_ret, aes(x = distance, y = bias, color = factor(subjectID))) +
  geom_point() +
  geom_line(aes(group = subjectID)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Groups ~ Region) +
  theme_minimal() +
  labs(color = "Subject ID") +
  ggtitle("Bias per Subject across Distance faceted by Group and Region (Retardance)")

# Figure 16: Bias per Subject Across Distance Faceted by Group and Region (Scattering)
bias_df <- means_mscat_data %>%
  mutate(
    pred_mean = fitted(bmod_comb_scat),       # posterior mean per row
    bias      = mean_value - pred_mean
  )

# Plot bias vs distance
ggplot(bias_df, aes(x = distance, y = bias, color = factor(subjectID))) +
  geom_point() +
  geom_line(aes(group = subjectID)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Groups ~ Region) +
  theme_minimal() +
  labs(color = "Subject ID") +
  ggtitle("Bias per Subject across Distance faceted by Group and Region (Scattering Coefficient)")


## adding random slopes for Region ##
# retardance
bmod_comb_ret_rs <- stan_glmer(
  mean_value ~ Groups * Region * distance_f + (1 + Region | subjectID),
  data = means_mret_data,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  cores = 4
)

# Figure 17
pp_check(bmod_comb_ret_rs)

# scattering
bmod_comb_scat_rs <- stan_glmer(
  mean_value ~ Groups * Region * distance_f +
    (1 + Region | subjectID),
  data = means_mscat_data,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  cores = 4
)

# Figure 18
pp_check(bmod_comb_scat_rs)

# Figure 19
plot(bmod_comb_ret_rs)

## Analysis for Retardance ##

post_ret <- as.data.frame(as.matrix(bmod_comb_ret_rs))
names(post_ret)

# FRONTAL CONTROL #
# identify all distance coef 
dist_terms <- grep("^distance_f", names(post_ret), value = TRUE)
dist_terms

# rename to specify frontal control
for (d in dist_terms) {
  new_name <- paste0("frontal_control_", d)
  post_ret[[new_name]] <- post_ret[[d]]
}

# rename
effects_cols <- grep("^frontal_control_distance_f", names(post_ret), value = TRUE)

summary_table_fc <- post_ret %>%
  select(all_of(effects_cols)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

print(summary_table_fc)

# FRONTAL EXPERIMENTAL #
for (d in dist_terms) {
  new_name <- paste0("frontal_experimental_", d)
  group_term <- paste0("Groupsexperimental:", d)   
  post_ret[[new_name]] <- post_ret[[d]] + post_ret[[group_term]]
}

effects_cols_fe <- grep("^frontal_experimental_distance_f", names(post_ret), value = TRUE)

# 4. Summarize posterior: mean + 95% CI
summary_table_fe <- post_ret %>%
  select(all_of(effects_cols_fe)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_fe

# OCCIPTIAL CONTROL #
for (d in dist_terms) {
  new_name <- paste0("occipital_control_", d)
  region_term <- paste0("Regionoccip:", d)  
  post_ret[[new_name]] <- post_ret[[d]] + post_ret[[region_term]]
}

effects_cols_oc <- grep("^occipital_control_distance_f", names(post_ret), value = TRUE)

summary_table_oc <- post_ret %>%
  select(all_of(effects_cols_oc)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_oc

# OCCIPITAL EXPERIMENTAL # 
for (d in dist_terms) {
  new_name <- paste0("occipital_experimental_", d)
  
  group_term  <- paste0("Groupsexperimental:", d)       
  region_term <- paste0("Regionoccip:", d)             
  inter_term  <- paste0("Groupsexperimental:Regionoccip:", d)  
  
  post_ret[[new_name]] <- post_ret[[d]] + 
    post_ret[[group_term]] + 
    post_ret[[region_term]] + 
    post_ret[[inter_term]]
}

effects_cols_oe <- grep("^occipital_experimental_distance_f", names(post_ret), value = TRUE)

summary_table_oe <- post_ret %>%
  select(all_of(effects_cols_oe)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_oe

# COMPARE CONTROL AND EXPERIMENTAL IN FRONTAL REGION # 
dist_terms <- grep("^Groupsexperimental:distance_f", names(post_ret), value = TRUE)

summary_diff <- data.frame(
  distance = gsub("Groupsexperimental:distance_f", "", dist_terms),
  mean = NA,
  ci_lower = NA,
  ci_upper = NA
)

for(i in seq_along(dist_terms)) {
  d <- dist_terms[i]
  posterior <- post_ret[[d]]
  
  summary_diff$mean[i] <- mean(posterior)
  summary_diff$ci_lower[i] <- quantile(posterior, 0.025)
  summary_diff$ci_upper[i] <- quantile(posterior, 0.975)
}

summary_diff

#COMPARE CONTROL AND EXPERIMENTAL IN OCCIPITAL REGION # 
dist_terms <- grep("^distance_f", names(post_ret), value = TRUE)

summary_diff_occip <- data.frame(
  distance = gsub("distance_f", "", dist_terms),
  mean = NA,
  ci_lower = NA,
  ci_upper = NA
)

for(d in dist_terms) {
  group_term <- paste0("Groupsexperimental:", d)
  interaction_term <- paste0("Groupsexperimental:Regionoccip:", d)
  
  posterior_diff <- post_ret[[group_term]] + post_ret[[interaction_term]]
  
  summary_diff_occip$mean[summary_diff_occip$distance == gsub("distance_f", "", d)] <- mean(posterior_diff)
  summary_diff_occip$ci_lower[summary_diff_occip$distance == gsub("distance_f", "", d)] <- quantile(posterior_diff, 0.025)
  summary_diff_occip$ci_upper[summary_diff_occip$distance == gsub("distance_f", "", d)] <- quantile(posterior_diff, 0.975)
}

summary_diff_occip




## Analysis for Scattering Coefficient ##

plot(bmod_comb_scat_rs)
post_scat <- as.data.frame(bmod_comb_scat_rs)

# FRONTAL CONTROL #
dist_terms_scat <- grep("^distance_f", names(post_scat), value = TRUE)
dist_terms_scat

for (d in dist_terms_scat) {
  new_name <- paste0("frontal_control_", d)
  post_scat[[new_name]] <- post_scat[[d]]
}

effects_scat_cols <- grep("^frontal_control_distance_f", names(post_scat), value = TRUE)

summary_table_fc_scat <- post_scat %>%
  select(all_of(effects_scat_cols)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_fc_scat

# FRONTAL EXPERIMENTAL #
for (d in dist_terms_scat) {
  new_name <- paste0("frontal_experimental_", d)
  group_term <- paste0("Groupsexperimental:", d)   
  post_scat[[new_name]] <- post_scat[[d]] + post_scat[[group_term]]
}

effects_scat_cols_fe <- grep("^frontal_experimental_distance_f", names(post_scat), value = TRUE)

summary_table_fe_scat <- post_scat %>%
  select(all_of(effects_scat_cols_fe)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_fe_scat

# OCCIPTIAL CONTROL #
for (d in dist_terms_scat) {
  new_name <- paste0("occipital_control_", d)
  region_term <- paste0("Regionoccip:", d)  
  post_scat[[new_name]] <- post_scat[[d]] + post_scat[[region_term]]
}

effects_scat_cols_oc <- grep("^occipital_control_distance_f", names(post_scat), value = TRUE)

summary_table_oc_scat <- post_scat %>%
  select(all_of(effects_scat_cols_oc)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_oc_scat

# OCCIPITAL EXPERIMENTAL # 
for (d in dist_terms_scat) {
  new_name <- paste0("occipital_experimental_", d)
  
  group_term  <- paste0("Groupsexperimental:", d)       
  region_term <- paste0("Regionoccip:", d)             
  inter_term  <- paste0("Groupsexperimental:Regionoccip:", d)  
  
  post_scat[[new_name]] <- post_scat[[d]] + 
    post_scat[[group_term]] + 
    post_scat[[region_term]] + 
    post_scat[[inter_term]]
}

effects_scat_cols_oe <- grep("^occipital_experimental_distance_f", names(post_scat), value = TRUE)

summary_table_oe_scat <- post_scat %>%
  select(all_of(effects_scat_cols_oe)) %>%
  summarise(across(everything(),
                   list(
                     mean = ~mean(.),
                     ci_lower = ~quantile(., 0.025),
                     ci_upper = ~quantile(., 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = c("effect", ".value"),
               names_pattern = "(.*)_(mean|ci_lower|ci_upper)")

summary_table_oe_scat

# COMPARE CONTROL AND EXPERIMENTAL IN FRONTAL REGION # 
summary_diff_scat <- data.frame(
  distance = gsub("Groupsexperimental:distance_f", "", dist_terms_scat),
  mean = NA,
  ci_lower = NA,
  ci_upper = NA
)

for(i in seq_along(dist_terms_scat)) {
  d <- dist_terms_scat[i]
  posterior <- post_scat[[d]]
  
  summary_diff_scat$mean[i] <- mean(posterior)
  summary_diff_scat$ci_lower[i] <- quantile(posterior, 0.025)
  summary_diff_scat$ci_upper[i] <- quantile(posterior, 0.975)
}

summary_diff_scat

# COMPARE CONTROL AND EXPERIMENTAL IN OCCIPITAL REGION # 
dist_terms <- grep("^distance_f", names(post_scat), value = TRUE)

summary_diff_occip_scat <- data.frame(
  distance = gsub("distance_f", "", dist_terms_scat),
  mean = NA,
  ci_lower = NA,
  ci_upper = NA
)

for(d in dist_terms_scat) {
  group_term <- paste0("Groupsexperimental:", d)
  interaction_term <- paste0("Groupsexperimental:Regionoccip:", d)
  
  posterior_diff <- post_scat[[group_term]] + post_scat[[interaction_term]]
  
  summary_diff_occip_scat$mean[summary_diff_occip_scat$distance == gsub("distance_f", "", d)] <- mean(posterior_diff)
  summary_diff_occip_scat$ci_lower[summary_diff_occip_scat$distance == gsub("distance_f", "", d)] <- quantile(posterior_diff, 0.025)
  summary_diff_occip_scat$ci_upper[summary_diff_occip_scat$distance == gsub("distance_f", "", d)] <- quantile(posterior_diff, 0.975)
}

summary_diff_occip_scat



