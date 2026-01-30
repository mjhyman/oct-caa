# load necessary libraries
library(nlme)
library(readxl)
library(rstanarm)
library(tidyverse)
library(tidybayes)
library(dplyr)
library(tidyr)
library(lme4)
library(shinystan)
library(broom)
library(tibble)
library(stringr)
library(ggplot2)


######## Output Filepath #######################################################
# Common directory
common_dir <- "/projectnb/npbssmic/ns/CAA"
input_filename <- "caa_all_radii_40um_donut_14-01-2026.xlsx"
input_filename <- "caa_all_radii_median_subtracted_40um_donut_15-01-2026.xlsx"
input_filename <- "caa_all_radii_percentage_diff_40um_donut_13-01-2026.xlsx"
input_file <- file.path(common_dir, input_filename)
output_filepath <- "/projectnb/npbssmic/ns/CAA/beta_stats/"
# Create a base filename
base_name <- tools::file_path_sans_ext(input_filename)
# Output filenames
mus_output_filename <- paste0(base_name,"_mus_summary_stats.csv")
mus_output_filename <- file.path(output_filepath, mus_output_filename)
ret_output_filename <- paste0(base_name,"_ret_summary_stats.csv")
ret_output_filename <- file.path(output_filepath, ret_output_filename)

######## 1) IMPORT DATA#########################################################
# retardance data
retardance_data <- read_excel(input_file, sheet = "retardance") %>% rename(retardance = OpticalProperty)
# scattering data
scattering_data <- read_excel(input_file, sheet = "scattering") %>% rename(scattering = OpticalProperty)

######## Create vector of mean at each distance ################################
# Compute per-group × distance raw summary:
# For each Groups (experimental/control), distance, and Region produce:
# - obs_mean: mean of measurements at that grouping
# - obs_sd: sample sd among measurements at that grouping
# - n: sample size used to compute the mean
# Also ensure distance is numeric

# Helper to normalize Region labels
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "Occip",
         ifelse(str_detect(r2, "front|frontal"), "Front", as.character(r)))
}

# Normalize region columns in raw data
retardance_raw <- retardance_data %>% mutate(Region = normalize_region(Region))
scattering_raw <- scattering_data %>% mutate(Region = normalize_region(Region))

raw_ret_summary <- retardance_raw %>%
  mutate(distance = as.numeric(distance)) %>%
  group_by(Region, Groups, distance) %>%
  summarise(
    obs_mean = mean(retardance, na.rm = TRUE),
    obs_sd = sd(retardance, na.rm = TRUE),
    n = sum(!is.na(retardance)),
    .groups = "drop"
  ) %>%
  rename(region = Region, group = Groups) %>%
  mutate(group = as.character(group))

raw_scat_summary <- scattering_raw %>%
  mutate(distance = as.numeric(distance)) %>%
  group_by(Region, Groups, distance) %>%
  summarise(
    obs_mean = mean(scattering, na.rm = TRUE),
    obs_sd = sd(scattering, na.rm = TRUE),
    n = sum(!is.na(scattering)),
    .groups = "drop"
  ) %>%
  rename(region = Region, group = Groups) %>%
  mutate(group = as.character(group))

### Now split into the four requested data frames.
# Retardance - Front
raw_ret_front <- raw_ret_summary %>%
  filter(region == "Front") %>%
  arrange(distance, group)
# Retardance - Occip
raw_ret_occip <- raw_ret_summary %>%
  filter(region == "Occip") %>%
  arrange(distance, group)
# Scattering - Front
raw_scat_front <- raw_scat_summary %>%
  filter(region == "Front") %>%
  arrange(distance, group)
# Scattering - Occip
raw_scat_occip <- raw_scat_summary %>%
  filter(region == "Occip") %>%
  arrange(distance, group)

# Quick sanity prints
message("Raw ret front rows: ", nrow(raw_ret_front))
message("Raw ret occip rows:  ", nrow(raw_ret_occip))
message("Raw scat front rows: ", nrow(raw_scat_front))
message("Raw scat occip rows:  ", nrow(raw_scat_occip))

# Example head() to inspect
head(raw_ret_front)
head(raw_ret_occip)
head(raw_scat_front)
head(raw_scat_occip)


######## Organize Retardance + Scattering Data  ################################

## average retardance across Groups, Region, subjectID,and distance
means_mret_data <- retardance_data %>%
  group_by(Groups, Stage, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(retardance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

## average scattering coefficient across Groups, Region, subjectID, and distance
means_mscat_data <- scattering_data %>%
  group_by(Groups, Stage, Region, subjectID, distance) %>%
  summarise(
    mean_value = mean(scattering, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distance_f = factor(distance))

# Ensure correct distances and factor levels
means_mret_data$distance <- factor(means_mret_data$distance, levels = seq(40, 480, 40))
means_mscat_data$distance <- factor(means_mscat_data$distance, levels = seq(40, 480, 40))

#### Split datasets into Front / Occip ####
# Retardance
means_mret_front <- filter(means_mret_data, Region == "front")
means_mret_occip  <- filter(means_mret_data, Region == "occip")
# Scattering
means_mscat_front <- filter(means_mscat_data, Region == "front")
means_mscat_occip  <- filter(means_mscat_data, Region == "occip")

######## Bayesian GLME w/ random slopes for Region + subID #####################
# Define the GLME model, excluding regions
glme_no_region <- as.formula("mean_value ~ Groups * Stage * distance_f + (1 | subjectID)")
glme_ran_grp <- as.formula("mean_value ~ Groups * Stage * distance_f + (1 + Groups | subjectID)")

# Stan GLME Sampling Settings
fit_ctrl <- list(chains = 4, iter = 4000, cores = 4, seed = 42)

### Fit stan_glmer models separately
# Retardance - Front
bmod_ret_front <- stan_glmer(
  glme_ran_grp,
  data = means_mret_front,
  family = gaussian(),
  chains = fit_ctrl$chains,
  iter = fit_ctrl$iter,
  cores = fit_ctrl$cores,
  seed = fit_ctrl$seed
)
# Retardance - Occip
bmod_ret_occip <- stan_glmer(
  glme_ran_grp,
  data = means_mret_occip,
  family = gaussian(),
  chains = fit_ctrl$chains,
  iter = fit_ctrl$iter,
  cores = fit_ctrl$cores,
  seed = fit_ctrl$seed
)
# Scattering - Front
bmod_scat_front <- stan_glmer(
  glme_ran_grp,
  data = means_mscat_front,
  family = gaussian(),
  chains = fit_ctrl$chains,
  iter = fit_ctrl$iter,
  cores = fit_ctrl$cores,
  seed = fit_ctrl$seed
)
# Scattering - Occip
bmod_scat_occip <- stan_glmer(
  glme_ran_grp,
  data = means_mscat_occip,
  family = gaussian(),
  chains = fit_ctrl$chains,
  iter = fit_ctrl$iter,
  cores = fit_ctrl$cores,
  seed = fit_ctrl$seed
)

# 4) Quick posterior checks (optional)
pp_check(bmod_ret_front);
pp_check(bmod_ret_occip)
pp_check(bmod_scat_front);
pp_check(bmod_scat_occip)



######### Posterior mean difference by distance (separate regions) #############
# Calculate posterior (Exp - Control) differences for each distance
# model: stan_glmer model (rstanarm)
# experimental_label: name of the experimental group coefficient, default "Groupsexperimental"
# region: NULL (ignore Region), or "Front"/"Occip" to include Region effects (for combined models)
calculate_posterior_diff <- function(model, experimental_label = "Groupsexperimental", region = NULL) {
  distances <- c(40, seq(80, 480, 40))
  
  # Extract draws matrix and convert to tibble with .draw index
  draws_mat <- as.matrix(model)
  draws_df <- as_tibble(draws_mat) %>% mutate(.draw = row_number())
  coef_names <- colnames(draws_mat)
  
  # helper to get coefficient vector or zeros if not present
  get_coef <- function(name) {
    if (!is.na(name) && name %in% colnames(draws_df)) draws_df[[name]] else rep(0, nrow(draws_df))
  }
  
  # Detect region coefficient name if present
  region_occ_name <- if ("Regionoccip" %in% coef_names) "Regionoccip" else NA_character_
  
  # Map distance main terms (distance_f40 is reference -> NA_character_)
  distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0("distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  
  # Map group x distance interaction names
  group_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0(experimental_label, ":distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  
  # Map region x distance names and group x region x distance names (if region used)
  region_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0("Regionoccip:distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  group_region_name <- if (paste0(experimental_label, ":Regionoccip") %in% coef_names) paste0(experimental_label, ":Regionoccip") else NA_character_
  group_region_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0(experimental_label, ":Regionoccip:distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  
  results <- list()
  preds_list <- list()
  
  for (d in distances) {
    # Build population-level prediction for control group at this distance
    pred_ctrl <- get_coef("(Intercept)")
    
    # Add distance main effect (absent for 40)
    if (!is.na(distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(distance_names[as.character(d)])
    
    # If region is specified and Occip, add region main and region:distance
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(region_occ_name)) pred_ctrl <- pred_ctrl + get_coef(region_occ_name)
      if (!is.na(region_distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(region_distance_names[as.character(d)])
    }
    
    # Build experimental prediction by adding group main and group interactions
    pred_exp <- pred_ctrl
    if (experimental_label %in% coef_names) pred_exp <- pred_exp + get_coef(experimental_label)
    if (!is.na(group_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(group_region_name)) pred_exp <- pred_exp + get_coef(group_region_name)
      if (!is.na(group_region_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_region_distance_names[as.character(d)])
    }
    
    # difference per draw
    diff_vec <- pred_exp - pred_ctrl
    
    preds_list[[as.character(d)]] <- tibble(
      distance = d,
      .draw = draws_df$.draw,
      diff = diff_vec
    )
  }
  
  all_diffs <- bind_rows(preds_list)
  
  summary_results <- all_diffs %>%
    group_by(distance) %>%
    summarise(
      posterior_mean = mean(diff, na.rm = TRUE),
      posterior_sd = sd(diff, na.rm = TRUE),
      n_draws = sum(!is.na(diff)),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(distance))
  
  return(list(summary = summary_results, raw = all_diffs))
}

### Calculate the posterior differences
# Retardance
res_ret_front <- calculate_posterior_diff(bmod_ret_front,
                                          experimental_label = "Groupsexperimental",
                                          region = NULL)
res_ret_occip <- calculate_posterior_diff(bmod_ret_occip,
                                          experimental_label = "Groupsexperimental",
                                          region = NULL)
# Scattering
res_scat_front <- calculate_posterior_diff(bmod_scat_front,
                                          experimental_label = "Groupsexperimental",
                                          region = NULL)
res_scat_occip <- calculate_posterior_diff(bmod_scat_occip,
                                          experimental_label = "Groupsexperimental",
                                          region = NULL)



#### Model predictions (mean + S.D.) at each distance  ####
# model: stan_glmer model
# experimental_label: name of experimental coefficient, e.g., "Groupsexperimental"
# region: NULL (no region), or "Front"/"Occip" to include Region terms when
#         present in the model
calculate_group_stats_by_distance <-
  function(model,
           experimental_label = "Groupsexperimental",
           region = NULL) {
  distances <- c(40, seq(80, 480, 40))
  
  draws_mat <- as.matrix(model)
  draws_df <- as_tibble(draws_mat) %>% mutate(.draw = row_number())
  coef_names <- colnames(draws_mat)
  
  # helper to return coefficient vector or zeros
  get_coef <- function(name) {
    if (!is.na(name) && name %in% colnames(draws_df)) draws_df[[name]] else rep(0, nrow(draws_df))
  }
  
  # region main name (if present)
  region_occ_name <- if ("Regionoccip" %in% coef_names) "Regionoccip" else NA_character_
  
  # distance main names (distance_f40 is reference -> NA)
  distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0("distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  
  # group x distance interaction names
  group_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0(experimental_label, ":distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  
  # region x distance and group x region interactions (if combined model)
  region_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0("Regionoccip:distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  group_region_name <- if (paste0(experimental_label, ":Regionoccip") %in% coef_names) paste0(experimental_label, ":Regionoccip") else NA_character_
  group_region_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0(experimental_label, ":Regionoccip:distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  
  # build raw per-draw predictions
  preds_list <- list()
  for (d in distances) {
    # control prediction (population-level)
    pred_ctrl <- get_coef("(Intercept)")
    if (!is.na(distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(region_occ_name)) pred_ctrl <- pred_ctrl + get_coef(region_occ_name)
      if (!is.na(region_distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(region_distance_names[as.character(d)])
    }
    
    # experimental prediction
    pred_exp <- pred_ctrl
    if (experimental_label %in% coef_names) pred_exp <- pred_exp + get_coef(experimental_label)
    if (!is.na(group_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(group_region_name)) pred_exp <- pred_exp + get_coef(group_region_name)
      if (!is.na(group_region_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_region_distance_names[as.character(d)])
    }
    
    preds_list[[as.character(d)]] <- tibble(
      distance = d,
      .draw = draws_df$.draw,
      pred_ctrl = pred_ctrl,
      pred_exp = pred_exp,
      diff = pred_exp - pred_ctrl
    )
  }
  
  raw_preds <- bind_rows(preds_list)
  
  # Summarize per distance: group means and sds and diff stats
  summary_results <- raw_preds %>%
    group_by(distance) %>%
    summarise(
      ctrl_mean = mean(pred_ctrl, na.rm = TRUE),
      ctrl_sd   = sd(pred_ctrl, na.rm = TRUE),
      exp_mean  = mean(pred_exp, na.rm = TRUE),
      exp_sd    = sd(pred_exp, na.rm = TRUE),
      diff_mean = mean(diff, na.rm = TRUE),
      diff_sd   = sd(diff, na.rm = TRUE),
      n_draws   = sum(!is.na(diff)),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(distance))
  
  return(list(summary = summary_results, raw = raw_preds))
}

## Retardance
# Frontal
res_ret_front_groupstats <-
  calculate_group_stats_by_distance(bmod_ret_front,
                                    experimental_label = "Groupsexperimental",
                                    region = NULL)
# Occipital
res_ret_occip_groupstats <-
  calculate_group_stats_by_distance(bmod_ret_occip,
                                    experimental_label = "Groupsexperimental",
                                    region = NULL)
print(res_ret_front_groupstats$summary, n = Inf)
print(res_ret_occip_groupstats$summary, n = Inf)

## Scattering
# Frontal
res_scat_front_groupstats <-
  calculate_group_stats_by_distance(bmod_scat_front,
                                    experimental_label = "Groupsexperimental",
                                    region = NULL)
# Occipital
res_scat_occip_groupstats <-
  calculate_group_stats_by_distance(bmod_scat_occip,
                                    experimental_label = "Groupsexperimental",
                                    region = NULL)
print(res_scat_front_groupstats$summary, n = Inf)
print(res_scat_occip_groupstats$summary, n = Inf)

#### Summarize raw diffs into mean, sd, and 95% credible interval ####
summarize_posterior_diffs <- function(raw_diffs_df, prob = 0.95) {
  alpha_low <- (1 - prob) / 2
  alpha_high <- 1 - alpha_low
  
  raw_diffs_df %>%
    group_by(region = if ("region" %in% colnames(.)) region else NA, distance) %>%
    summarise(
      posterior_mean = mean(diff, na.rm = TRUE),
      posterior_sd = sd(diff, na.rm = TRUE),
      lower = quantile(diff, probs = alpha_low, na.rm = TRUE),
      upper = quantile(diff, probs = alpha_high, na.rm = TRUE),
      n_draws = sum(!is.na(diff)),
      .groups = "drop"
    ) %>%
    # if region column is NA (not present), set to "All"
    mutate(region = ifelse(is.na(region), "All", as.character(region))) %>%
    arrange(region, as.numeric(distance))
}

### Print Summary Stats and Plot mean +/- 95% CI across distances
# summary_df: result from summarize_posterior_diffs
# facet_by_region: TRUE -> separate panels per region; FALSE -> overlay regions
# fixed_y: if TRUE uses same y-scale across facets; if FALSE uses free y-scales
plot_diffs_with_ci <- function(summary_df, title = "Posterior mean difference (Exp - Control)",
                               facet_by_region = FALSE, fixed_y = TRUE,
                               point_size = 1.5, line_size = 1) {
  # ensure distance numeric
  summary_df <- summary_df %>% mutate(distance = as.numeric(distance))
  
  p <- ggplot(summary_df, aes(x = distance, y = posterior_mean, color = region, fill = region)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(aes(group = region), size = line_size) +
    geom_point(size = point_size) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_continuous(breaks = seq(40, 480, 40)) +
    labs(x = "Distance", y = "Posterior mean (Exp - Control)", title = title, color = "Region", fill = "Region") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (facet_by_region) {
    scales_arg <- if (fixed_y) "fixed" else "free_y"
    p <- p + facet_wrap(~ region, nrow = 1, scales = scales_arg) +
      theme(legend.position = "none")  # legend redundant in faceted view
  }
  
  return(p)
}

### Retardance - Frontal
# Printed Results
summary_ret_front <- summarize_posterior_diffs(res_ret_front$raw)
print(''); print('Results for Frontal Retardance')
print(summary_ret_front, n = Inf)
# Graphical Results
p_ret_front <- plot_diffs_with_ci(summary_ret_front,
                                  title = "Retardance (Front): Exp - Control",
                                  facet_by_region = FALSE)
print(p_ret_front)

### Retardance - Occipital
# Printed Results
summary_ret_occip <- summarize_posterior_diffs(res_ret_occip$raw)
print(''); print('Results for Occipital Retardance')
print(summary_ret_occip, n = Inf)
# Graphical Results
p_ret_occip <- plot_diffs_with_ci(summary_ret_occip,
                                  title = "Retardance (Occip): Exp - Control",
                                  facet_by_region = FALSE)
print(p_ret_occip)

### Scattering - Frontal
# Printed Results
summary_scat_front <- summarize_posterior_diffs(res_scat_front$raw)
print(''); print('Results for Frontal Scattering')
print(summary_scat_front, n = Inf)
# Graphical Results
p_scat_front <- plot_diffs_with_ci(summary_scat_front,
                                  title = "Scattering (Front): Exp - Control",
                                  facet_by_region = FALSE)
print(p_scat_front)

### Scattering - Occipital
# Printed Results
summary_scat_occip <- summarize_posterior_diffs(res_scat_occip$raw)
print(''); print('Results for Occipital Scattering')
print(summary_scat_occip, n = Inf)
# Graphical Results
p_scat_occip <- plot_diffs_with_ci(summary_scat_occip,
                                  title = "Scattering (Occip): Exp - Control",
                                  facet_by_region = FALSE)
print(p_scat_occip)

### Save all summaries to CSV


######## Posterior Mean Differences (Exp - Control) #########################
# # Set the group labels as variables
# experimental_group <- "Groupsexperimental"  # Coefficient for the experimental group
# # No need to directly reference "(Intercept)" when manually handling random effects
# 
# ### Define function to calculate posterior mean difference for a model
# calculate_posterior_diff <- function(model, experimental_label) {
#   group_samples <- model %>%
#     # Correctly access the intercept
#     spread_draws(!!sym(experimental_label), `(Intercept)`) %>%
#     
#     # Calculate the mean difference: Experimental - Control
#     # Use `(Intercept)` for direct access
#     mutate(mean_diff = !!sym(experimental_label) - `(Intercept)`) %>%
#     select(.draw, mean_diff)
#   
#   return(group_samples)
# }


######### Posterior mean differnce by region + distance ##########
# ### Scattering
# # Calculate posterior mean differences for the scattering model
# posterior_scat_diff <- calculate_posterior_diff(bmod_comb_scat_rs,
#                                                 experimental_group)
# # Summarize the results
# summary_scat <- posterior_scat_diff %>%
#   summarise(mean_diff = mean(mean_diff),
#             sd_diff = sd(mean_diff))
# # Print summary
# print(summary_scat)
# 
# ### Retardance
# # Calculate posterior mean differences for the scattering model
# posterior_ret_diff <- calculate_posterior_diff(bmod_comb_ret_rs,
#                                                 experimental_group)
# # Summarize the results
# summary_ret <- posterior_ret_diff %>%
#   summarise(mean_diff = mean(mean_diff),
#             sd_diff = sd(mean_diff))
# 
# # Print summary
# print(summary_ret)

#### Function to calculate posterior mean differences by region and distance####
# calculate_posterior_diff_by_region_distance <- function(model, experimental_label = "Groupsexperimental") {
#   distances <- c(40, seq(80, 480, 40))
#   regions <- c("Front", "Occip")
#   
#   draws_mat <- as.matrix(model)
#   draws_df <- as_tibble(draws_mat) %>% mutate(.draw = row_number())
#   coef_names <- colnames(draws_mat)
#   
#   # Helper to pull coef vector or zeros if missing
#   get_coef <- function(name) {
#     if (!is.na(name) && name %in% colnames(draws_df)) draws_df[[name]] else rep(0, nrow(draws_df))
#   }
#   
#   intercept <- "(Intercept)"
#   region_occ <- if ("Regionoccip" %in% coef_names) "Regionoccip" else NA_character_
#   
#   # Map distance names (distance_f40 absent intentionally)
#   distance_names <- setNames(
#     sapply(distances, function(d) {
#       nm <- paste0("distance_f", d)
#       if (nm %in% coef_names) nm else NA_character_
#     }, USE.NAMES = FALSE),
#     distances
#   )
#   
#   # Group x distance names
#   group_distance_names <- setNames(
#     sapply(distances, function(d) {
#       nm <- paste0(experimental_label, ":distance_f", d)
#       if (nm %in% coef_names) nm else NA_character_
#     }, USE.NAMES = FALSE),
#     distances
#   )
#   
#   # Group x Region and Group x Region x distance names
#   group_region <- if (paste0(experimental_label, ":Regionoccip") %in% coef_names) paste0(experimental_label, ":Regionoccip") else NA_character_
#   group_region_distance_names <- setNames(
#     sapply(distances, function(d) {
#       nm <- paste0(experimental_label, ":Regionoccip:distance_f", d)
#       if (nm %in% coef_names) nm else NA_character_
#     }, USE.NAMES = FALSE),
#     distances
#   )
#   
#   # Region x distance names
#   region_distance_names <- setNames(
#     sapply(distances, function(d) {
#       nm <- paste0("Regionoccip:distance_f", d)
#       if (nm %in% coef_names) nm else NA_character_
#     }, USE.NAMES = FALSE),
#     distances
#   )
#   
#   results_list <- list()
#   for (reg in regions) {
#     for (d in distances) {
#       # control prediction
#       pred_ctrl <- get_coef(intercept)
#       if (!is.na(distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(distance_names[as.character(d)])
#       if (reg == "Occip") {
#         if (!is.na(region_occ)) pred_ctrl <- pred_ctrl + get_coef(region_occ)
#         if (!is.na(region_distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(region_distance_names[as.character(d)])
#       }
#       # experimental prediction
#       pred_exp <- pred_ctrl
#       if (experimental_label %in% coef_names) pred_exp <- pred_exp + get_coef(experimental_label)
#       if (!is.na(group_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_distance_names[as.character(d)])
#       if (reg == "Occip" && !is.na(group_region)) pred_exp <- pred_exp + get_coef(group_region)
#       if (reg == "Occip" && !is.na(group_region_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_region_distance_names[as.character(d)])
#       
#       diff_vec <- pred_exp - pred_ctrl
#       results_list[[paste0(reg, "_", d)]] <- tibble(region = reg, distance = d, .draw = draws_df$.draw, diff = diff_vec)
#     }
#   }
#   
#   all_diffs <- bind_rows(results_list)
#   summary_results <- all_diffs %>%
#     group_by(region, distance) %>%
#     summarise(posterior_mean = mean(diff, na.rm = TRUE),
#               posterior_sd = sd(diff, na.rm = TRUE),
#               n_draws = sum(!is.na(diff)),
#               .groups = "drop") %>%
#     arrange(region, distance)
#   list(summary = summary_results, raw = all_diffs)
# }

# ### Calculate posterior mean difference by region and distance
# mus_post_diff <- calculate_posterior_diff_by_region_distance(bmod_comb_scat_rs)
# ret_post_diff <- calculate_posterior_diff_by_region_distance(bmod_comb_ret_rs)
# ### Print scattering and retardance summaries
# print(mus_post_diff$summary, n=Inf)
# print(ret_post_diff$summary, n=Inf)

########## Summaries for combined region GLME ##################################
# ##### 1) summarize raw per-draw diffs (mean, sd, 95% credible interval)####
# summarize_posterior_diffs <- function(raw_diffs_df, prob = 0.95) {
#   alpha_low <- (1 - prob) / 2
#   alpha_high <- 1 - alpha_low
#   raw_diffs_df %>%
#     group_by(region, distance) %>%
#     summarise(
#       posterior_mean = mean(diff, na.rm = TRUE),
#       posterior_sd = sd(diff, na.rm = TRUE),
#       lower = quantile(diff, probs = alpha_low, na.rm = TRUE),
#       upper = quantile(diff, probs = alpha_high, na.rm = TRUE),
#       n_draws = sum(!is.na(diff)),
#       .groups = "drop"
#     ) %>%
#     arrange(region, distance)
# }
# 
# ############ 2) Plot posterior mean with ribbon for credible interval ##########
# plot_posterior_diffs_facet <- function(summary_df, title = "Posterior mean difference (Experimental - Control)", facet_by = "region", horizontal = TRUE) {
#   summary_df <- summary_df %>% mutate(distance = as.numeric(distance))
#   
#   p <- ggplot(summary_df, aes(x = distance, y = posterior_mean)) +
#     geom_line(aes(color = region), size = 1) +
#     geom_point(aes(color = region), size = 1.5) +
#     geom_ribbon(aes(ymin = lower, ymax = upper, fill = region), alpha = 0.2, color = NA) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
#     scale_x_continuous(breaks = seq(40, 480, 40)) +
#     labs(x = "Distance", y = "Posterior mean (Exp - Control)",
#          title = title, color = "Region", fill = "Region") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   # facet: horizontal (panels side-by-side) or vertical (stacked)
#   if (horizontal) {
#     p <- p + facet_wrap(~ region, nrow = 1, scales = "free_y")
#   } else {
#     p <- p + facet_wrap(~ region, ncol = 1, scales = "free_y")
#   }
#   
#   return(p)
# }
# 
# # # Create summaries with 95% intervals
# # summary_mus <- summarize_posterior_diffs(mus_post_diff$raw, prob = 0.95)
# # summary_ret  <- summarize_posterior_diffs(ret_post_diff$raw, prob = 0.95)
# # 
# # # Print the summaries
# # print(summary_mus, n = Inf)
# # print(summary_ret, n = Inf)
# # 
# # # Scattering Coefficient w/ vertical stacking:
# # p_mus_vertical <- plot_posterior_diffs_facet(summary_mus,
# #                                               title = "Scattering: Posterior mean difference",
# #                                               horizontal = FALSE)
# # print(p_mus_vertical)
# # 
# # # Retardance w/ vertical stacking:
# # p_ret_vertical <- plot_posterior_diffs_facet(summary_ret,
# #                                               title = "Retardance: Posterior mean difference",
# #                                               horizontal = FALSE)
# # print(p_ret_vertical)
# 
# ############### Save summaries to CSV ##########################################
# # Write CSV files
# # write.csv(summary_mus, mus_output_filename, row.names = FALSE)
# # write.csv(summary_ret, ret_output_filename, row.names = FALSE)
