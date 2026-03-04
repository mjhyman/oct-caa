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
library(posterior)
library(purrr)


# TODO: output the summary statistics for supplemental

######## Filepaths #######################################################
# Common directory
common_dir <- "/projectnb/npbssmic/ns/CAA"
input_filenames <- c(
  "caa_all_radii_40um_donut_14-01-2026.xlsx",
  "caa_all_radii_median_subtracted_40um_donut_15-01-2026.xlsx",
  "caa_all_radii_percentage_diff_40um_donut_13-01-2026.xlsx"
)

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

#### Helper to normalize Region labels ####
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "Occip",
         ifelse(str_detect(r2, "front|frontal"), "Front", as.character(r)))
}

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


# Function to create posterior summary table from a fitted model
# model: fitted model object (rstanarm, brms, or stanfit) where as.matrix(model) works
# probs: interval probabilities (defaults to 0.025 and 0.975 for 95% CI)
# use_median: if TRUE uses posterior median as point estimate; if FALSE uses posterior mean
# exclude_pars: optional character vector of parameter name regexes to exclude (e.g., "b\\[" for random effects)
# n_digits: digits for rounding in printed table
# save_csv: optional path to save CSV
# Required packages
posterior_summary_table <- function(model,
                                    probs = c(0.025, 0.975),
                                    use_median = TRUE,
                                    exclude_pars = c(), # regex patterns to exclude
                                    n_digits = 3,
                                    save_csv = NULL) {
  # Try to get draws with chain structure using posterior
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for reliable R-hat and ESS. Please install: install.packages('posterior')")
  }
  draws_arr <- tryCatch(posterior::as_draws_array(model),
                        error = function(e) {
                          stop("Failed to convert model to draws array via posterior::as_draws_array(): ", e$message)
                        })
  # draws_arr: iterations x chains x parameters
  # Get parameter names in the array
  param_names <- dimnames(draws_arr)$variable
  if (is.null(param_names)) param_names <- colnames(as.matrix(model))
  
  # Optionally exclude parameters by patterns
  if (!is.null(exclude_pars) && length(exclude_pars) > 0) {
    exclude_mask <- rep(FALSE, length(param_names))
    for (pat in exclude_pars) exclude_mask <- exclude_mask | stringr::str_detect(param_names, pat)
    keep_idx <- which(!exclude_mask)
    if (length(keep_idx) == 0) stop("No parameters left after applying exclude_pars")
    draws_arr <- draws_arr[, , keep_idx, drop = FALSE]
    param_names <- param_names[keep_idx]
  }
  
  # Use posterior::summarize_draws for efficient summaries (includes r_hat and ess_bulk)
  sums <- posterior::summarize_draws(draws_arr) # tibble with variable, mean, sd, ess_bulk, r_hat, quantiles...
  # Ensure sums rows correspond to param_names
  # summarize_draws returns rows in same order as variables in draws_arr
  # Extract needed columns; some older posterior versions might name ess differently; check column names
  has_rhat <- "rhat" %in% colnames(sums)
  has_ess <- "ess_bulk" %in% colnames(sums)
  
  # Extract quantiles according to probs
  qnames <- paste0("q", round(probs * 100, 0))
  # posterior::summarize_draws stores quantiles in columns named with probs, e.g., '2.5%' or 'q025' depending on version.
  # safer to compute quantiles directly from draws_arr per parameter.
  
  # Convert draws_arr to matrix per parameter to compute quantiles reliably:
  draws_mat <- posterior::as_draws_matrix(draws_arr) # draws x params, preserves order
  # Now build table
  summaries <- purrr::map_dfr(seq_along(param_names), function(i) {
    pname <- param_names[i]
    vec <- draws_mat[, i]
    vec <- vec[is.finite(vec)]
    est <- if (use_median) median(vec) else mean(vec)
    sd_x <- sd(vec)
    ci_low <- quantile(vec, probs[1], na.rm = TRUE)
    ci_high <- quantile(vec, probs[2], na.rm = TRUE)
    # find rhat and ess from sums (match by variable)
    rhat_val <- if (has_rhat && pname %in% sums$variable) {
      as.numeric(sums$rhat[which(sums$variable == pname)])
    } else NA_real_
    neff_val <- if (has_ess && pname %in% sums$variable) {
      as.numeric(sums$ess_bulk[which(sums$variable == pname)])
    } else NA_real_
    tibble(
      parameter = pname,
      estimate = est,
      sd = sd_x,
      ci_lower = ci_low,
      ci_upper = ci_high,
      n_eff = neff_val,
      Rhat = rhat_val
    )
  })
  
  # Round numeric columns for presentation
  num_cols <- c("estimate", "sd", "ci_lower", "ci_upper", "n_eff", "Rhat")
  summaries <- summaries %>%
    mutate(across(all_of(num_cols), ~ ifelse(is.na(.), NA, signif(., digits = n_digits))))
  
  # Save CSV if requested
  if (!is.null(save_csv)) {
    if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
    readr::write_csv(summaries, save_csv)
  }
  
  return(summaries)
}


#### Iterate over each spreadsheet
for (fname in input_filenames) {
  
  print('ANALYZING DATA FOR')
  print(fname)
  input_file <- file.path(common_dir, fname)

  ######## 1) IMPORT DATA#######################################################
  # retardance data
  retardance_data <- read_excel(input_file, sheet = "retardance") %>% rename(retardance = OpticalProperty)
  # scattering data
  scattering_data <- read_excel(input_file, sheet = "scattering") %>% rename(scattering = OpticalProperty)
  
  ######## Create vector of mean at each distance ##############################
  # Compute per-group × distance raw summary:
  # For each Groups (experimental/control), distance, and Region produce:
  # - obs_mean: mean of measurements at that grouping
  # - obs_sd: sample sd among measurements at that grouping
  # - n: sample size used to compute the mean
  # Also ensure distance is numeric
  
  # Normalize region columns in raw data
  retardance_raw <- retardance_data %>% mutate(Region = normalize_region(Region))
  scattering_raw <- scattering_data %>% mutate(Region = normalize_region(Region))
  
  # Summarize retardance data
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
  
  # Summarize scattering data
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
  
  
  ######## Organize Retardance + Scattering Data  ##############################
  
  ## average retardance across Groups, Region, subjectID,and distance
  means_mret_data <- retardance_data %>%
    group_by(Groups, Stage, Region, subjectID, distance) %>%
    summarise(
      mean_value = mean(retardance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(distance_f = factor(distance))
  
  ## average MUS across Groups, Region, subjectID, and distance
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
  
  
  ############ Calculate group statistics by distance ###################
  
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
  summary_mus_front <- summarize_posterior_diffs(res_scat_front$raw)
  print(''); print('Results for Frontal Scattering')
  print(summary_mus_front, n = Inf)
  # Graphical Results
  p_scat_front <- plot_diffs_with_ci(summary_mus_front,
                                    title = "Scattering (Front): Exp - Control",
                                    facet_by_region = FALSE)
  print(p_scat_front)
  
  ### Scattering - Occipital
  # Printed Results
  summary_mus_occip <- summarize_posterior_diffs(res_scat_occip$raw)
  print(''); print('Results for Occipital Scattering')
  print(summary_mus_occip, n = Inf)
  # Graphical Results
  p_scat_occip <- plot_diffs_with_ci(summary_mus_occip,
                                    title = "Scattering (Occip): Exp - Control",
                                    facet_by_region = FALSE)
  print(p_scat_occip)
  
  #### Save all summaries to CSV ####
  # Create full output folder name
  output_filepath <- "/projectnb/npbssmic/ns/CAA/beta_stats"
  # Create a base filename
  base_name <- tools::file_path_sans_ext(fname)
  # MUS Output filenames
  mus_front_fout <- paste0(base_name,"__mus_front_summary_stats.csv")
  mus_front_fout <- file.path(output_filepath, base_name, mus_front_fout)
  mus_occip_fout <- paste0(base_name,"__mus_occip_summary_stats.csv")
  mus_occip_fout <- file.path(output_filepath, base_name, mus_occip_fout)
  # RET Output filenames
  ret_front_fout <- paste0(base_name,"__ret_front_summary_stats.csv")
  ret_front_fout <- file.path(output_filepath, base_name, ret_front_fout)
  ret_occip_fout <- paste0(base_name,"__ret_occip_summary_stats.csv")
  ret_occip_fout <- file.path(output_filepath, base_name, ret_occip_fout)
  
  # Save output files
  write.csv(res_scat_front_groupstats$summary, file=mus_front_fout, row.names = FALSE)
  write.csv(res_scat_occip_groupstats$summary, file=mus_occip_fout, row.names = FALSE)
  write.csv(res_ret_front_groupstats$summary, file=ret_front_fout, row.names = FALSE)
  write.csv(res_ret_occip_groupstats$summary, file=ret_occip_fout, row.names = FALSE)
  
  ### Posteriors
  # Scattering Frontal
  mus_front_fout <- paste0(base_name,"__mus_front_posterior.csv")
  mus_front_fout <- file.path(output_filepath, base_name, mus_front_fout)
  # Scattering Occipital
  mus_occip_fout <- paste0(base_name,"__mus_occip_posterior.csv")
  mus_occip_fout <- file.path(output_filepath, base_name, mus_occip_fout)
  # Retardance Frontal
  ret_front_fout <- paste0(base_name,"__ret_front_posterior.csv")
  ret_front_fout <- file.path(output_filepath, base_name, ret_front_fout)
  # Retardance Occipital
  ret_occip_fout <- paste0(base_name,"__ret_occip_posterior.csv")
  ret_occip_fout <- file.path(output_filepath, base_name, ret_occip_fout)
  # Posterior for Frontal Scattering
  tbl <- posterior_summary_table(
    bmod_scat_front,
    probs = c(0.025, 0.975),
    use_median = TRUE,
    exclude_pars = c("^b\\[", "^Sigma", "^Sigma\\[", "^b\\("), # exclude per-subject random-effect draws etc.
    n_digits = 3,
    save_csv = mus_front_fout
  )
  # Posterior for Occipital Scattering
  tbl <- posterior_summary_table(
    bmod_scat_occip,
    probs = c(0.025, 0.975),
    use_median = TRUE,
    exclude_pars = c("^b\\[", "^Sigma", "^Sigma\\[", "^b\\("), # exclude per-subject random-effect draws etc.
    n_digits = 3,
    save_csv = mus_occip_fout
  )
  # Posterior for Frontal Retardance
  tbl <- posterior_summary_table(
    bmod_ret_front,
    probs = c(0.025, 0.975),
    use_median = TRUE,
    exclude_pars = c("^b\\[", "^Sigma", "^Sigma\\[", "^b\\("), # exclude per-subject random-effect draws etc.
    n_digits = 3,
    save_csv = ret_front_fout
  )
  # Posterior for Occipital Retardance
  tbl <- posterior_summary_table(
    bmod_ret_occip,
    probs = c(0.025, 0.975),
    use_median = TRUE,
    exclude_pars = c("^b\\[", "^Sigma", "^Sigma\\[", "^b\\("), # exclude per-subject random-effect draws etc.
    n_digits = 3,
    save_csv = ret_occip_fout
  )

}
