# load necessary libraries
library(readxl)
library(rstanarm)
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(posterior)
library(purrr)
library(readr)


######## Filepaths ###########################################################
common_dir <- "/projectnb/npbssmic/ns/CAA"
input_filenames <- c(
  "caa_all_radii_percentage_diff_40um_donut_5-Mar-2026.xlsx"
)


######### Helper: normalize Region labels ####################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}


######### Helper: fit one within-stage model #################################
# Within a fixed optical property x region x source x stage, the only predictor
# is distance. Subjects per stage: 1 for control/mild/moderate, 2 for severe.
# Because no stage has >= 3 subjects we always use stan_glm. When 2 subjects are
# present, subjectID is added as a fixed covariate; the fitted curve is later
# averaged over the two subjects.
fit_stage_model <- function(data, response, fit_ctrl) {
  n_subj <- length(unique(data$subjectID))
  f <- if (n_subj >= 2) {
    as.formula(paste0(response, " ~ distance_f + subjectID"))
  } else {
    as.formula(paste0(response, " ~ distance_f"))
  }
  stan_glm(
    f, data = data, family = gaussian(),
    chains = fit_ctrl$chains, iter = fit_ctrl$iter,
    cores  = fit_ctrl$cores,  seed = fit_ctrl$seed
  )
}


######### Helper: posterior fitted means by distance #########################
# Returns a draws x distance matrix of the model-implied mean at each distance,
# averaged over whatever subjects are present in the fitting data. Uses only the
# distances actually present in this stage's data (so posterior_epred never has
# to predict a factor level that was not estimated).
posterior_fitted_by_distance <- function(model, data) {
  d_present    <- sort(unique(data$distance))
  subj_present <- unique(data$subjectID)
  
  newdata <- expand.grid(
    distance  = d_present,
    subjectID = subj_present,
    stringsAsFactors = FALSE
  )
  newdata$distance_f <- factor(newdata$distance, levels = seq(40, 480, 40))
  
  ep <- posterior_epred(model, newdata = newdata)   # draws x nrow(newdata)
  
  # average over subject columns within each distance -> draws x ndist
  draws_by_d <- sapply(d_present, function(dd) {
    cols <- which(newdata$distance == dd)
    rowMeans(ep[, cols, drop = FALSE])
  })
  if (is.null(dim(draws_by_d))) draws_by_d <- matrix(draws_by_d, nrow = nrow(ep))
  colnames(draws_by_d) <- as.character(d_present)
  
  list(draws = draws_by_d, distances = d_present)
}


######### Method 1: model-based difference (stage k - stage 0) ###############
# Differences the two independent stage models' posterior fitted curves draw by
# draw, then summarizes the difference distribution at each common distance.
model_diff_table <- function(model0, data0, modelk, datak, prob = 0.95) {
  f0 <- posterior_fitted_by_distance(model0, data0)
  fk <- posterior_fitted_by_distance(modelk, datak)
  
  common <- sort(intersect(f0$distances, fk$distances))
  if (length(common) == 0) return(NULL)
  
  nd      <- min(nrow(f0$draws), nrow(fk$draws))  # align draw counts (normally equal)
  a_low   <- (1 - prob) / 2
  a_high  <- 1 - a_low
  
  purrr::map_dfr(common, function(dd) {
    v0   <- f0$draws[seq_len(nd), as.character(dd)]
    vk   <- fk$draws[seq_len(nd), as.character(dd)]
    diff <- vk - v0
    tibble(
      distance  = dd,
      diff_mean = mean(diff),
      diff_sd   = sd(diff),
      lower     = unname(quantile(diff, a_low)),
      upper     = unname(quantile(diff, a_high)),
      n_draws   = length(diff)
    )
  }) %>% arrange(distance)
}


######### Method 2: raw-mean difference (stage k - stage 0) ##################
# Difference of raw per-distance means, with the propagated standard error of
# the difference of means: SE = sqrt(s_k^2/n_k + s_0^2/n_0). The raw components
# (s, n) for each stage are also returned so a pooled SD sqrt(s_k^2 + s_0^2) can
# be recomputed if that convention is preferred.
raw_means_by_distance <- function(data, response) {
  data %>%
    group_by(distance) %>%
    summarise(
      mean = mean(.data[[response]], na.rm = TRUE),
      sd   = sd(.data[[response]],   na.rm = TRUE),
      n    = sum(!is.na(.data[[response]])),
      .groups = "drop"
    )
}

raw_diff_table <- function(data0, datak, response) {
  r0 <- raw_means_by_distance(data0, response) %>% rename(mean_0 = mean, sd_0 = sd, n_0 = n)
  rk <- raw_means_by_distance(datak, response) %>% rename(mean_k = mean, sd_k = sd, n_k = n)
  
  inner_join(r0, rk, by = "distance") %>%
    mutate(
      diff    = mean_k - mean_0,
      se_diff = sqrt(sd_k^2 / n_k + sd_0^2 / n_0)
    ) %>%
    arrange(distance) %>%
    select(distance, mean_0, mean_k, diff, sd_0, sd_k, n_0, n_k, se_diff)
}


######### Helper: overlay plot of both methods ###############################
plot_stage_diffs <- function(model_df, raw_df, title) {
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_ribbon(data = model_df,
                aes(x = distance, ymin = lower, ymax = upper),
                alpha = 0.2, fill = "steelblue") +
    geom_line(data = model_df,
              aes(x = distance, y = diff_mean), color = "steelblue", linewidth = 0.8) +
    geom_point(data = model_df,
               aes(x = distance, y = diff_mean), color = "steelblue", size = 1.3) +
    geom_errorbar(data = raw_df,
                  aes(x = distance, ymin = diff - se_diff, ymax = diff + se_diff),
                  color = "firebrick", width = 8, alpha = 0.7) +
    geom_point(data = raw_df,
               aes(x = distance, y = diff), color = "firebrick", size = 1.3) +
    facet_wrap(~ stage_label, nrow = 1) +
    scale_x_continuous(breaks = seq(40, 480, 80)) +
    labs(
      x = "Distance", y = "Difference (Stage k - Stage 0)",
      title = title,
      subtitle = "blue = model fitted-curve difference (95% CI ribbon); red = raw-mean difference (\u00b1 SE)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8))
}


######## BEGIN LOGIC #########################################################

fname      <- input_filenames[1]
input_file <- file.path(common_dir, fname)
message("ANALYZING DATA FOR: ", fname)

######## 1) IMPORT DATA ######################################################
# Groups encodes the measurement source:
#   experimental = measurement around the enlarged perivascular space (EPVS)
#   control      = measurement around the vessel
# Stage encodes CAA severity: 0 = control, 1 = mild, 2 = moderate, 3 = severe.
retardance_data <- read_excel(input_file, sheet = "retardance") %>%
  rename(retardance = OpticalProperty) %>%
  mutate(
    Region = normalize_region(Region),
    Stage  = as.integer(as.numeric(Stage)),
    Groups = factor(Groups, levels = c("control", "experimental"))
  )

scattering_data <- read_excel(input_file, sheet = "scattering") %>%
  rename(scattering = OpticalProperty) %>%
  mutate(
    Region = normalize_region(Region),
    Stage  = as.integer(as.numeric(Stage)),
    Groups = factor(Groups, levels = c("control", "experimental"))
  )

######## 2) DEFINE LOOP PARAMETERS ###########################################
properties     <- list(retardance = retardance_data, scattering = scattering_data)
regions        <- c("front", "occip")
sources        <- c(epvs = "experimental", vessel = "control")   # name -> Groups level
nonzero_stages <- c(1L, 2L, 3L)
stage_labels   <- c("0" = "control", "1" = "mild", "2" = "moderate", "3" = "severe")

fit_ctrl        <- list(chains = 4, iter = 4000, cores = 4, seed = 42)
output_filepath <- "/projectnb/npbssmic/ns/CAA/beta_stats"
base_name       <- tools::file_path_sans_ext(fname)

######## 3) MAIN LOOP: PROPERTY x REGION x SOURCE x (STAGE k - STAGE 0) ######
for (prop in names(properties)) {
  pdata <- properties[[prop]]
  
  for (rgn in regions) {
    for (src in names(sources)) {
      grp <- sources[[src]]
      
      message("\n===== ", prop, " | ", rgn, " | ", src, " (Groups = ", grp, ") =====")
      
      base_sub <- pdata %>%
        filter(Region == rgn, Groups == grp) %>%
        mutate(distance_f = factor(distance, levels = seq(40, 480, 40)))
      
      data0 <- base_sub %>% filter(Stage == 0)
      if (nrow(data0) == 0) {
        message("  Skipping – no stage 0 (control) data for this property/region/source")
        next
      }
      
      # Fit the stage 0 model once and reuse for every nonzero-stage comparison
      message("  Fitting stage 0 (control) model...")
      model0 <- fit_stage_model(data0, prop, fit_ctrl)
      
      out_dir <- file.path(output_filepath, base_name, "stage_vs_control", rgn, src)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      model_accum <- list()
      raw_accum   <- list()
      
      for (k in nonzero_stages) {
        datak <- base_sub %>% filter(Stage == k)
        if (nrow(datak) == 0) {
          message("  Skipping stage ", k, " (", stage_labels[as.character(k)],
                  ") – no data for this property/region/source")
          next
        }
        
        message("  Fitting stage ", k, " (", stage_labels[as.character(k)], ") model...")
        modelk <- fit_stage_model(datak, prop, fit_ctrl)
        
        ## Method 1: model fitted-curve difference
        mdiff <- model_diff_table(model0, data0, modelk, datak)
        ## Method 2: raw-mean difference
        rdiff <- raw_diff_table(data0, datak, prop)
        
        if (is.null(mdiff) || nrow(rdiff) == 0) {
          message("  No common distances between stage ", k, " and stage 0 – skipping")
          next
        }
        
        ## Save both summary tables
        model_fout <- file.path(
          out_dir,
          paste0(base_name, "__", prop, "_", rgn, "_", src,
                 "_stage", k, "_vs_stage0_model_diff.csv")
        )
        raw_fout <- file.path(
          out_dir,
          paste0(base_name, "__", prop, "_", rgn, "_", src,
                 "_stage", k, "_vs_stage0_raw_diff.csv")
        )
        write.csv(mdiff, model_fout, row.names = FALSE)
        write.csv(rdiff, raw_fout,   row.names = FALSE)
        
        ## Print summary stats
        message("  --- Method 1 (model) : stage ", k, " - stage 0 ---")
        print(mdiff, n = Inf)
        message("  --- Method 2 (raw)   : stage ", k, " - stage 0 ---")
        print(rdiff, n = Inf)
        
        ## Accumulate for the overlay plot
        sl <- paste0("Stage ", k, " (", stage_labels[as.character(k)], ") - control")
        model_accum[[as.character(k)]] <- mdiff %>% mutate(stage = k, stage_label = sl)
        raw_accum[[as.character(k)]]   <- rdiff %>% mutate(stage = k, stage_label = sl)
      }
      
      ## Overlay plot across all available nonzero stages
      if (length(model_accum) > 0) {
        model_all <- bind_rows(model_accum)
        raw_all   <- bind_rows(raw_accum)
        p <- plot_stage_diffs(
          model_all, raw_all,
          title = paste0(prop, " (", rgn, ", ", src, "): stage - control")
        )
        print(p)
      }
      
      message("  Finished: ", out_dir)
    }
  }
}