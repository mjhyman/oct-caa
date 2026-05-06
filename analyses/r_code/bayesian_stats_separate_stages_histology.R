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

######## Configuration & Path Setup ############################################
common_dir <- "/projectnb/npbssmic/ns/CAA/histology/stats/"
input_filename <- "histogram_matched_donuts_2026-May-05.xlsx"
input_file <- file.path(common_dir, input_filename)

# Updated discrete distances
dist_levels <- c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)
stain_tabs <- c("lhe", "gfap", "cd68")

# Generate Timestamped Directory
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_root <- file.path(common_dir, paste0("bayesian_stats_histo_", timestamp))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Response variable name update
response_var <- "stain_value"

######### Helper: normalize Region labels #####################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}

######### Helper: formula selection ##########################################
make_formula <- function(response, data) {
  n_subjects <- length(unique(data$subjectID))
  if (n_subjects >= 3) {
    as.formula(paste0(response, " ~ Groups * distance_f + (1 + Groups | subjectID)"))
  } else {
    # Fallback for small N
    as.formula(paste0(response, " ~ Groups * distance_f"))
  }
}

######### Calculation Helpers #################################################
calculate_posterior_diff <- function(model, experimental_label = "GroupsExperimental") {
  # Uses global dist_levels
  draws_mat  <- as.matrix(model)
  draws_df   <- as_tibble(draws_mat) %>% mutate(.draw = row_number())
  coef_names <- colnames(draws_mat)
  
  get_coef <- function(name) {
    if (!is.na(name) && name %in% colnames(draws_df)) draws_df[[name]] else rep(0, nrow(draws_df))
  }
  
  preds_list <- list()
  for (d in dist_levels) {
    d_nm <- paste0("distance_f", d)
    int_nm <- paste0(experimental_label, ":distance_f", d)
    
    # Simple contrast: Effect of Group at Distance D
    diff_vec <- get_coef(experimental_label) + get_coef(int_nm)
    
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
      lower = quantile(diff, 0.025),
      upper = quantile(diff, 0.975),
      .groups = "drop"
    )
  
  return(list(summary = summary_results, raw = all_diffs))
}

######## Initialize Readme #####################################################
# Summarizes analysis details [cite: 2, 42, 56]
readme_path <- file.path(output_root, "readme.txt")
writeLines(paste0("Analysis Run: ", timestamp, 
                  "\nInput: ", input_filename, 
                  "\nDistances: ", paste(dist_levels, collapse=", "),
                  "\nStains processed: ", paste(stain_tabs, collapse=", "),
                  "\nNotes: Values represent averaged ROI within histogram-matched/z-score normalized sections.",
                  "\n------------------------------------------\n"), readme_path)

######## MAIN EXECUTION LOOP ###################################################
fit_ctrl <- list(chains = 4, iter = 4000, cores = 4, seed = 42)

for (stain in stain_tabs) {
  message("\nProcessing Stain Tab: ", stain)
  
  # Import specific tab with updated column names and levels [cite: 32, 44, 46, 52, 54]
  raw_data <- tryCatch({
    read_excel(input_file, sheet = stain) %>%
      mutate(
        Region = normalize_region(Region),
        Stage  = as.integer(as.numeric(Stage)),
        # Match new capitalization
        Groups = factor(Groups, levels = c("Control", "Experimental")),
        # Match new capitalized Distance column
        distance_f = factor(Distance, levels = dist_levels)
      )
  }, error = function(e) {
    message("Could not read sheet: ", stain)
    return(NULL)
  })
  
  if (is.null(raw_data)) next
  
  stages  <- c(0L, 1L, 2L, 3L)
  regions <- c("front", "occip")
  
  for (stg in stages) {
    for (rgn in regions) {
      
      sub_data <- raw_data %>% filter(Stage == stg, Region == rgn)
      
      if (nrow(sub_data) < 5) next # Skip if insufficient data
      
      label <- paste0(stain, "_stage", stg, "_", rgn)
      stain_dir <- file.path(output_root, stain, paste0("stage", stg))
      dir.create(stain_dir, recursive = TRUE, showWarnings = FALSE)
      
      message("  Fitting model for ", label)
      
      # Use updated response variable "stain_value" [cite: 45, 51]
      formula <- make_formula(response_var, sub_data)
      n_subs <- length(unique(sub_data$subjectID))
      
      bmod <- if(n_subs >= 3) {
        stan_glmer(formula, data = sub_data, chains=fit_ctrl$chains, iter=fit_ctrl$iter, seed=fit_ctrl$seed, refresh=0)
      } else {
        stan_glm(formula, data = sub_data, chains=fit_ctrl$chains, iter=fit_ctrl$iter, seed=fit_ctrl$seed, refresh=0)
      }
      
      # Calculate stats and differences [cite: 27, 31]
      res <- calculate_posterior_diff(bmod, experimental_label = "GroupsExperimental")
      
      # Save Summary CSV
      write.csv(res$summary, file.path(stain_dir, paste0(label, "_diff_summary.csv")), row.names = FALSE)
      
      # Save Visualization
      p <- ggplot(res$summary, aes(x = as.numeric(as.character(distance)), y = posterior_mean)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "darkred") +
        geom_line(color = "darkred", size = 1) +
        geom_point() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(title = paste(stain, "| Stage", stg, "|", rgn),
             subtitle = "Posterior Difference (Experimental - Control)",
             x = "Distance (microns)", y = "Mean Difference (stain_value)") +
        theme_minimal()
      
      ggsave(file.path(stain_dir, paste0(label, "_plot.png")), p, width = 8, height = 5)
      
      # Append to Readme [cite: 42, 56]
      cat(paste0("Processed: ", label, " | N_Subjects: ", n_subs, " | Rows: ", nrow(sub_data), "\n"), 
          file = readme_path, append = TRUE)
    }
  }
}

message("\nAnalysis Complete. Results saved to: ", output_root)