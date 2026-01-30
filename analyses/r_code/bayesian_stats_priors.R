# load necessary libraries
library(nlme)
library(readxl)
library(rstanarm)
library(tidybayes)
library(tidyverse)
library(lme4)
library(loo)
library(shinystan)


######## IMPORT DATA############################################################
### RETARDANCE ###
# retardance data set
retardance_data <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_percentage_diff_40um_donut_13-01-2026.xlsx", sheet = "retardance")
# rename OpticalProperty to retardance
retardance_data <- retardance_data %>% rename(retardance = OpticalProperty)
# Create a copy
m_retardance_data <- retardance_data

### SCATTERING ###
# import data
scattering_data <- read_excel("/projectnb/npbssmic/ns/CAA/caa_all_radii_percentage_diff_40um_donut_13-01-2026.xlsx", sheet = "scattering")
# rename OpticalProperty to be scattering 
scattering_data <- scattering_data %>% rename(scattering = OpticalProperty)
# Create a copy
m_scattering_data <- scattering_data

######## Organize Retardance + Scattering Data  ################################

## Combine Retardance into a single table (front + occip) ##
# average retardance across Groups, Region, subjectID,and distance
means_mret_data <- m_retardance_data %>%
  group_by(Groups, Stage, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(retardance, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

## Combine Scattering into a single table (front + occip) ##
# average scattering coefficient across Groups, Region, subjectID, and distance
means_mscat_data <- m_scattering_data %>%
  group_by(Groups, Stage, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(scattering, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

######## Type Conversions + Standardize ########################################
prep_data <- function(df) {
  df %>%
    mutate(
      group = factor(Groups),              # experimental vs control
      region = factor(Region),             # frontal / occipital
      subjectID = factor(subjectID),
      stage_num = as.integer(Stage),       # 0..3 ordinal numeric
      # standardize stage and distance for modeling
      stage_c = (stage_num - mean(stage_num, na.rm = TRUE)) / sd(stage_num, na.rm = TRUE),
      distance_num = as.numeric(distance),
      distance_c = (distance_num - mean(distance_num, na.rm = TRUE)) / sd(distance_num, na.rm = TRUE),
      # numeric 0/1 contrast for region (explicit for random slope)
      region_num = ifelse(region == levels(region)[1], 0L, 1L)
    )
}
# Apply function to dataset
means_mret_data <- prep_data(means_mret_data)
means_mscat_data <- prep_data(means_mscat_data)

######## Outcome summaries for priors ##########################################
sd_ret <- sd(means_mret_data$mean_value, na.rm = TRUE)
median_ret <- median(means_mret_data$mean_value, na.rm = TRUE)
cat("Retardance: median =", median_ret, " sd =", sd_ret, "\n")

sd_scat <- sd(means_mscat_data$mean_value, na.rm = TRUE)
median_scat <- median(means_mscat_data$mean_value, na.rm = TRUE)
cat("Scattering: median =", median_scat, " sd =", sd_scat, "\n")

######## Create scaled response ################################################
means_mret_data <- means_mret_data %>%
  mutate(mean_value_sc = (mean_value - mean(mean_value, na.rm = TRUE)) / sd_ret)
means_mscat_data <- means_mscat_data %>%
  mutate(mean_value_sc = (mean_value - mean(mean_value, na.rm = TRUE)) / sd_scat)


######## Prior definitions (data-informed) #####################################
# Fixed effects prior
fixed_scale_ret <- 0.5 * sd_ret
fixed_scale_scat <- 0.5 * sd_scat

# Random effects prior (exp w/ mean = 0.5 * sd_outcome -> rate = 1 / mean)
re_rate_ret <- 1 / (0.5 * sd_ret)
re_rate_scat <- 1 / (0.5 * sd_scat)

# Residual prior (exponential with mean = sd_outcome -> rate = 1 / sd_outcome)
res_rate_ret <- 1 / sd_ret
res_rate_scat <- 1 / sd_scat

######## Model definitions #####################################################
# Full 4-way interaction as original (but stage and distance are numeric scaled)
formula_full_raw <- mean_value ~ group * stage_c * region * distance_c +
                    (1 + region_num | subjectID)
formula_full_sc  <- mean_value_sc ~ group * stage_c * region * distance_c +
                    (1 + region_num | subjectID)
# Alternative: stage as ordered factor
formula_stage_factor <- mean_value ~ group * factor(stage_num, ordered = TRUE) *
                    region * distance_c + (1 + region_num | subjectID)

######## Fit Models (raw scale w/ data-informed priors)  #######################
#Retardance raw-scale model
bmod_comb_ret_rs <- stan_glmer(
  formula = formula_full_raw,
  data = means_mret_data,
  family = gaussian(),
  prior = normal(0, fixed_scale_ret, autoscale = FALSE),    # fixed effects
  prior_intercept = normal(median_ret, sd_ret, autoscale = FALSE),
  prior_covariance = decov(regularization = 2,
                           concentration = sd_ret,
                           prior_sds = exponential(rate = re_rate_ret)),
  prior_aux = exponential(rate = res_rate_ret),
  chains = 4,
  iter = 3000,
  cores = 4,
  seed = 42,
  refresh = 200
)

#Scattering raw-scale model
bmod_comb_scat_rs <- stan_glmer(
  formula = formula_full_raw,
  data = means_mscat_data,
  family = gaussian(),
  prior = normal(0, fixed_scale_scat, autoscale = FALSE),
  prior_intercept = normal(median_scat, sd_scat, autoscale = FALSE),
  prior_covariance = decov(regularization = 2,
                           concentration = sd_scat,
                           prior_sds = exponential(rate = re_rate_scat)),
  prior_aux = exponential(rate = res_rate_scat),
  chains = 4,
  iter = 3000,
  cores = 4,
  seed = 42,
  refresh = 200
)

######## Fit Models (standardized response scale)  #############################
# Retardance
bmod_comb_ret_sc <- stan_glmer(
  formula = formula_full_sc,
  data = means_mret_data,
  family = gaussian(),
  prior = normal(0, 0.5, autoscale = FALSE),
  prior_intercept = normal(0, 1, autoscale = FALSE),
  prior_covariance = decov(regularization = 2,
                           concentration = 1,
                           prior_sds = exponential(rate = 1 / 0.5)),
  prior_aux = exponential(rate = 1),
  chains = 4,
  iter = 3000,
  cores = 4,
  seed = 42,
  refresh = 200
)

# Scattering
bmod_comb_scat_sc <- stan_glmer(
  formula = formula_full_sc,
  data = means_mscat_data,
  family = gaussian(),
  prior = normal(0, 0.5, autoscale = FALSE),
  prior_intercept = normal(0, 1, autoscale = FALSE),
  prior_covariance = decov(regularization = 2,
                           concentration = 1,
                           prior_sds = exponential(rate = 1 / 0.5)),
  prior_aux = exponential(rate = 1),
  chains = 4,
  iter = 3000,
  cores = 4,
  seed = 42,
  refresh = 200
)

######## Fit Alternative Model #################################################
bmod_stage_factor_ret <- stan_glmer(
  formula = formula_stage_factor,
  data = means_mret_data,
  family = gaussian(),
  prior = normal(0, fixed_scale_ret, autoscale = FALSE),
  prior_intercept = normal(median_ret, sd_ret, autoscale = FALSE),
  prior_covariance = decov(regularization = 2,
                           concentration = sd_ret,
                           prior_sds = exponential(rate = re_rate_ret)),
  prior_aux = exponential(rate = res_rate_ret),
  chains = 4,
  iter = 3000,
  cores = 4,
  seed = 42,
  refresh = 200
)

#### TODO: add rest from saving, diagnostics, etc.


######## Bayesian GLME w/ random slopes for Region + subID #####################

# Set the random seed
set.seed(42)

# Retardance
bmod_comb_ret_rs <- stan_glmer(
  mean_value ~ Groups * Stage * Region * distance_f + (1 + Region | subjectID),
  data = means_mret_data,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  cores = 4
)

# Figure 17 - RETARDANCE posterior comparison
pp_check(bmod_comb_ret_rs)

# scattering
bmod_comb_scat_rs <- stan_glmer(
  mean_value ~ Groups * Stage * Region * distance_f + (1 + Region | subjectID),
  data = means_mscat_data,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  cores = 4
)

# Figure 18 - MUS posterior comparison
pp_check(bmod_comb_scat_rs)

# Figure 19
plot(bmod_comb_ret_rs)

######## RETARDANCE SUMMARY TABLES #############################################

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



######## SCATTERING SUMMARY TABLES##############################################
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

# Print to console
summary_diff_occip_scat
