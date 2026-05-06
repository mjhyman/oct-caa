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
library(readr)


######## Filepaths #######################################################
common_dir <- "/projectnb/npbssmic/ns/CAA/histology/stats"
input_filenames <- c(
  "histogram_matched_donuts_2026-May-05.xlsx"
)
# Generate Timestamped Directory [cite: 1, 39, 55]
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_root <- file.path(common_dir, paste0("bayesian_stats_histo_", timestamp))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
distance <- c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)

######### Helper: normalize Region labels #####################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}


######### Helper: choose formula for per-stage models ########################
make_formula <- function(response, data) {
  n_subjects <- length(unique(data$subjectID))
  if (n_subjects >= 3) {
    as.formula(paste0(response, " ~ Groups * distance_f + (1 + Groups | subjectID)"))
  } else if (n_subjects == 2) {
    message("  2 subjects found – fitting fixed-effects model with subjectID covariate")
    as.formula(paste0(response, " ~ Groups * distance_f * subjectID"))
  } else {
    message("  1 subject found – fitting fixed-effects-only model")
    as.formula(paste0(response, " ~ Groups * distance_f"))
  }
}


######### Helper: choose formula for combined (all-stage) models #############
make_formula_combined <- function(response, data) {
  n_subjects <- length(unique(data$subjectID))
  if (n_subjects >= 3) {
    as.formula(paste0(response, " ~ Groups * distance_f + (1 + distance_f | subjectID)"))
  } else if (n_subjects == 2) {
    message("  2 subjects found – fitting fixed-effects model with subjectID covariate")
    as.formula(paste0(response, " ~ Groups * distance_f * subjectID"))
  } else {
    message("  1 subject found – fitting fixed-effects-only model")
    as.formula(paste0(response, " ~ Groups * distance_f"))
  }
}


######### Posterior mean difference by distance ###############################
calculate_posterior_diff <- function(model, experimental_label = "GroupsExperimental", region = NULL) {
  distances <- c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)
  
  draws_mat  <- as.matrix(model)
  draws_df   <- as_tibble(draws_mat) %>% mutate(.draw = row_number())
  coef_names <- colnames(draws_mat)
  
  get_coef <- function(name) {
    if (!is.na(name) && name %in% colnames(draws_df)) draws_df[[name]] else rep(0, nrow(draws_df))
  }
  
  region_occ_name <- if ("Regionoccip" %in% coef_names) "Regionoccip" else NA_character_
  
  distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0("distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  group_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0(experimental_label, ":distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
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
  
  preds_list <- list()
  for (d in distances) {
    pred_ctrl <- get_coef("(Intercept)")
    if (!is.na(distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(region_occ_name)) pred_ctrl <- pred_ctrl + get_coef(region_occ_name)
      if (!is.na(region_distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(region_distance_names[as.character(d)])
    }
    pred_exp <- pred_ctrl
    if (experimental_label %in% coef_names) pred_exp <- pred_exp + get_coef(experimental_label)
    if (!is.na(group_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(group_region_name)) pred_exp <- pred_exp + get_coef(group_region_name)
      if (!is.na(group_region_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_region_distance_names[as.character(d)])
    }
    diff_vec <- pred_exp - pred_ctrl
    preds_list[[as.character(d)]] <- tibble(
      distance = d,
      .draw    = draws_df$.draw,
      diff     = diff_vec
    )
  }
  
  all_diffs <- bind_rows(preds_list)
  
  summary_results <- all_diffs %>%
    group_by(distance) %>%
    summarise(
      posterior_mean = mean(diff, na.rm = TRUE),
      posterior_sd   = sd(diff,   na.rm = TRUE),
      n_draws        = sum(!is.na(diff)),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(distance))
  
  return(list(summary = summary_results, raw = all_diffs))
}


######### Model predictions (mean + SD) at each distance #####################
calculate_group_stats_by_distance <- function(model, experimental_label = "GroupsExperimental", region = NULL) {
  distances <- c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)
  
  draws_mat  <- as.matrix(model)
  draws_df   <- as_tibble(draws_mat) %>% mutate(.draw = row_number())
  coef_names <- colnames(draws_mat)
  
  get_coef <- function(name) {
    if (!is.na(name) && name %in% colnames(draws_df)) draws_df[[name]] else rep(0, nrow(draws_df))
  }
  
  region_occ_name <- if ("Regionoccip" %in% coef_names) "Regionoccip" else NA_character_
  
  distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0("distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
  group_distance_names <- setNames(
    sapply(distances, function(d) {
      nm <- paste0(experimental_label, ":distance_f", d)
      if (nm %in% coef_names) nm else NA_character_
    }, USE.NAMES = FALSE),
    distances
  )
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
  
  preds_list <- list()
  for (d in distances) {
    pred_ctrl <- get_coef("(Intercept)")
    if (!is.na(distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(region_occ_name)) pred_ctrl <- pred_ctrl + get_coef(region_occ_name)
      if (!is.na(region_distance_names[as.character(d)])) pred_ctrl <- pred_ctrl + get_coef(region_distance_names[as.character(d)])
    }
    pred_exp <- pred_ctrl
    if (experimental_label %in% coef_names) pred_exp <- pred_exp + get_coef(experimental_label)
    if (!is.na(group_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_distance_names[as.character(d)])
    if (!is.null(region) && tolower(region) == "occip") {
      if (!is.na(group_region_name)) pred_exp <- pred_exp + get_coef(group_region_name)
      if (!is.na(group_region_distance_names[as.character(d)])) pred_exp <- pred_exp + get_coef(group_region_distance_names[as.character(d)])
    }
    preds_list[[as.character(d)]] <- tibble(
      distance  = d,
      .draw     = draws_df$.draw,
      pred_ctrl = pred_ctrl,
      pred_exp  = pred_exp,
      diff      = pred_exp - pred_ctrl
    )
  }
  
  raw_preds <- bind_rows(preds_list)
  
  summary_results <- raw_preds %>%
    group_by(distance) %>%
    summarise(
      ctrl_mean = mean(pred_ctrl, na.rm = TRUE),
      ctrl_sd   = sd(pred_ctrl,   na.rm = TRUE),
      exp_mean  = mean(pred_exp,  na.rm = TRUE),
      exp_sd    = sd(pred_exp,    na.rm = TRUE),
      diff_mean = mean(diff,      na.rm = TRUE),
      diff_sd   = sd(diff,        na.rm = TRUE),
      n_draws   = sum(!is.na(diff)),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(distance))
  
  return(list(summary = summary_results, raw = raw_preds))
}


######### Summarize raw diffs into mean, SD, and 95% credible interval ########
summarize_posterior_diffs <- function(raw_diffs_df, prob = 0.95) {
  alpha_low  <- (1 - prob) / 2
  alpha_high <- 1 - alpha_low
  
  raw_diffs_df %>%
    group_by(region = if ("region" %in% colnames(.)) region else NA, distance) %>%
    summarise(
      posterior_mean = mean(diff, na.rm = TRUE),
      posterior_sd   = sd(diff,   na.rm = TRUE),
      lower          = quantile(diff, probs = alpha_low,  na.rm = TRUE),
      upper          = quantile(diff, probs = alpha_high, na.rm = TRUE),
      n_draws        = sum(!is.na(diff)),
      .groups = "drop"
    ) %>%
    mutate(region = ifelse(is.na(region), "All", as.character(region))) %>%
    arrange(region, as.numeric(distance))
}


######### Plot mean +/- 95% CI across distances ###############################
plot_diffs_with_ci <- function(summary_df,
                               title = "Posterior mean difference (Exp - Control)",
                               facet_by_region = FALSE, fixed_y = TRUE,
                               point_size = 1.5, line_size = 1) {
  summary_df <- summary_df %>% mutate(distance = as.numeric(distance))
  
  p <- ggplot(summary_df, aes(x = distance, y = posterior_mean, color = region, fill = region)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(aes(group = region), size = line_size) +
    geom_point(size = point_size) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_continuous(breaks = c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)) +
    labs(x = "Distance", y = "Posterior mean (Exp - Control)", title = title,
         color = "Region", fill = "Region") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (facet_by_region) {
    scales_arg <- if (fixed_y) "fixed" else "free_y"
    p <- p + facet_wrap(~ region, nrow = 1, scales = scales_arg) +
      theme(legend.position = "none")
  }
  return(p)
}


######### Posterior summary table #############################################
posterior_summary_table <- function(model,
                                    probs        = c(0.025, 0.975),
                                    use_median   = TRUE,
                                    exclude_pars = c(),
                                    n_digits     = 3,
                                    save_csv     = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required.")
  
  draws_arr <- tryCatch(
    posterior::as_draws_array(model),
    error = function(e) stop("Failed to convert model: ", e$message)
  )
  param_names <- dimnames(draws_arr)$variable
  if (is.null(param_names)) param_names <- colnames(as.matrix(model))
  
  if (!is.null(exclude_pars) && length(exclude_pars) > 0) {
    exclude_mask <- rep(FALSE, length(param_names))
    for (pat in exclude_pars) exclude_mask <- exclude_mask | stringr::str_detect(param_names, pat)
    keep_idx <- which(!exclude_mask)
    if (length(keep_idx) == 0) stop("No parameters left after applying exclude_pars")
    draws_arr   <- draws_arr[, , keep_idx, drop = FALSE]
    param_names <- param_names[keep_idx]
  }
  
  sums      <- posterior::summarize_draws(draws_arr)
  has_rhat  <- "rhat"     %in% colnames(sums)
  has_ess   <- "ess_bulk" %in% colnames(sums)
  draws_mat <- posterior::as_draws_matrix(draws_arr)
  
  summaries <- purrr::map_dfr(seq_along(param_names), function(i) {
    pname    <- param_names[i]
    vec      <- draws_mat[, i]
    vec      <- vec[is.finite(vec)]
    est      <- if (use_median) median(vec) else mean(vec)
    ci_low   <- quantile(vec, probs[1], na.rm = TRUE)
    ci_hi    <- quantile(vec, probs[2], na.rm = TRUE)
    rhat_val <- if (has_rhat && pname %in% sums$variable) as.numeric(sums$rhat[sums$variable     == pname]) else NA_real_
    neff_val <- if (has_ess  && pname %in% sums$variable) as.numeric(sums$ess_bulk[sums$variable == pname]) else NA_real_
    tibble(parameter = pname, estimate = est, sd = sd(vec),
           ci_lower = ci_low, ci_upper = ci_hi, n_eff = neff_val, Rhat = rhat_val)
  })
  
  num_cols  <- c("estimate", "sd", "ci_lower", "ci_upper", "n_eff", "Rhat")
  summaries <- summaries %>%
    mutate(across(all_of(num_cols), ~ ifelse(is.na(.), NA, signif(., digits = n_digits))))
  
  if (!is.null(save_csv)) readr::write_csv(summaries, save_csv)
  return(summaries)
}


######### Helper: fit, summarize, and save one per-stage modality #############
fit_and_save_modality <- function(data, response, stg, rgn, label, out_dir,
                                  base_name, fit_ctrl) {
  if (nrow(data) == 0 || length(unique(data$Groups)) < 2) {
    message("  Skipping ", response, " – insufficient data for ", label)
    return(invisible(NULL))
  }
  
  message("  Fitting ", response, " model...")
  n_subjects <- length(unique(data$subjectID))
  formula    <- make_formula(response, data)
  
  # stan_glmer requires random effects terms; use stan_glm for <3 subjects
  if (n_subjects < 3) {
    bmod <- stan_glm(formula, data = data, family = gaussian(),
                     chains  = fit_ctrl$chains, 
                     iter    = fit_ctrl$iter,
                     warmup  = fit_ctrl$warmup,
                     cores   = fit_ctrl$cores,  
                     seed    = fit_ctrl$seed,
                     control = fit_ctrl$control)
  } else {
    bmod <- stan_glmer(formula, data = data, family = gaussian(),
                       chains  = fit_ctrl$chains, 
                       iter    = fit_ctrl$iter,
                       warmup  = fit_ctrl$warmup,
                       cores   = fit_ctrl$cores,  
                       seed    = fit_ctrl$seed,
                       control = fit_ctrl$control)
  }
  
  print(pp_check(bmod) + ggtitle(paste("PP check –", response, "–", label)))
  
  res    <- calculate_posterior_diff(bmod, experimental_label = "GroupsExperimental", region = NULL)
  gs     <- calculate_group_stats_by_distance(bmod, experimental_label = "GroupsExperimental", region = NULL)
  sum_df <- summarize_posterior_diffs(res$raw)
  
  message("  Results – ", response, " – ", label)
  print(sum_df, n = Inf)
  print(plot_diffs_with_ci(sum_df,
                           title = paste0(response, " (", rgn, ", Stage ", stg, "): Exp - Control")))
  
  write.csv(gs$summary,
            file.path(out_dir, paste0(base_name, "__", response, "_", label, "_summary_stats.csv")),
            row.names = FALSE)
  posterior_summary_table(
    bmod,
    probs        = c(0.025, 0.975),
    use_median   = TRUE,
    exclude_pars = c("^b\\[", "^Sigma", "^Sigma\\[", "^b\\("),
    n_digits     = 3,
    save_csv     = file.path(out_dir, paste0(base_name, "__", response, "_", label, "_posterior.csv"))
  )
  
  invisible(bmod)
}


######### Helper: fit, summarize, and save one combined modality ##############
fit_and_save_modality_combined <- function(data, response, rgn, out_dir,
                                           base_name, fit_ctrl) {
  label <- paste0("combined_", rgn)
  
  if (nrow(data) == 0 || length(unique(data$Groups)) < 2) {
    message("  Skipping ", response, " – insufficient data for ", label)
    return(invisible(NULL))
  }
  
  message("  Fitting ", response, " model...")
  n_subjects <- length(unique(data$subjectID))
  formula    <- make_formula_combined(response, data)
  
  if (n_subjects < 3) {
    bmod <- stan_glm(formula, data = data, family = gaussian(),
                     chains = fit_ctrl$chains, iter = fit_ctrl$iter,
                     cores  = fit_ctrl$cores,  seed = fit_ctrl$seed)
  } else {
    bmod <- stan_glmer(formula, data = data, family = gaussian(),
                       chains = fit_ctrl$chains, iter = fit_ctrl$iter,
                       cores  = fit_ctrl$cores,  seed = fit_ctrl$seed)
  }
  
  print(pp_check(bmod) + ggtitle(paste("PP check –", response, "–", label)))
  
  res    <- calculate_posterior_diff(bmod, experimental_label = "GroupsExperimental", region = NULL)
  gs     <- calculate_group_stats_by_distance(bmod, experimental_label = "GroupsExperimental", region = NULL)
  sum_df <- summarize_posterior_diffs(res$raw)
  
  message("  Results – ", response, " – ", label)
  print(sum_df, n = Inf)
  print(plot_diffs_with_ci(sum_df,
                           title = paste0(response, " (", rgn, ", All Stages): Exp - Control")))
  
  write.csv(gs$summary,
            file.path(out_dir, paste0(base_name, "__", response, "_", label, "_summary_stats.csv")),
            row.names = FALSE)
  posterior_summary_table(
    bmod,
    probs        = c(0.025, 0.975),
    use_median   = TRUE,
    exclude_pars = c("^b\\[", "^Sigma", "^Sigma\\[", "^b\\("),
    n_digits     = 3,
    save_csv     = file.path(out_dir, paste0(base_name, "__", response, "_", label, "_posterior.csv"))
  )
  
  invisible(bmod)
}


######## BEGIN LOGIC ###########################################################

fname      <- input_filenames[1]
input_file <- file.path(common_dir, fname)
message("ANALYZING DATA FOR: ", fname)

######## 1) IMPORT DATA ########################################################
lhe_data <- read_excel(input_file, sheet = "lhe") %>%
  rename(lhe = stain_value, distance = Distance) %>%
  mutate(
    Region = normalize_region(Region),
    Stage  = as.integer(as.numeric(Stage)),
    Groups = factor(Groups, levels = c("Control", "Experimental"))
  )

gfap_data <- read_excel(input_file, sheet = "gfap") %>%
  rename(gfap = stain_value, distance = Distance) %>%
  mutate(
    Region = normalize_region(Region),
    Stage  = as.integer(as.numeric(Stage)),
    Groups = factor(Groups, levels = c("Control", "Experimental"))
  )

cd68_data <- read_excel(input_file, sheet = "cd68") %>%
  rename(cd68 = stain_value, distance = Distance) %>%
  mutate(
    Region = normalize_region(Region),
    Stage  = as.integer(as.numeric(Stage)),
    Groups = factor(Groups, levels = c("Control", "Experimental"))
  )

######## 2) DEFINE LOOP PARAMETERS ############################################
stages  <- c(0L, 1L, 2L, 3L)
regions <- c("front", "occip")

fit_ctrl        <- list(
  chains = 4,
  iter = 10000,
  warmup = 5000,
  cores = 4,
  seed = 42,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
    )
  )
base_name       <- tools::file_path_sans_ext(fname)

######## 3) MAIN LOOP: PER STAGE × REGION #####################################
for (stg in stages) {
  for (rgn in regions) {
    
    label <- paste0("stage", stg, "_", rgn)
    message("\n===== Processing: ", label, " =====")
    
    lhe_sub  <- lhe_data %>%
      filter(Stage == stg, Region == rgn) %>%
      mutate(distance_f = factor(distance, levels = c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)))
    gfap_sub <- gfap_data %>%
      filter(Stage == stg, Region == rgn) %>%
      mutate(distance_f = factor(distance, levels = c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)))
    cd68_sub <- cd68_data %>%
      filter(Stage == stg, Region == rgn) %>%
      mutate(distance_f = factor(distance, levels = c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)))
    
    out_dir <- file.path(output_root, paste0("stage", stg), label)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    fit_and_save_modality(lhe_sub, "lhe", stg, rgn, label, out_dir, base_name, fit_ctrl)
    fit_and_save_modality(gfap_sub,"gfap",stg, rgn, label, out_dir, base_name, fit_ctrl)
    fit_and_save_modality(cd68_sub,"cd68",stg, rgn, label, out_dir, base_name, fit_ctrl)
    
    message("  Finished: ", out_dir)
  }
}