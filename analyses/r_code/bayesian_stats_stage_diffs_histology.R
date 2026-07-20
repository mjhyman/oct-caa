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


######## Filepaths ###########################################################
common_dir <- "/projectnb/npbssmic/ns/CAA/histology/stats"
input_filenames <- c(
  "histogram_matched_donuts_2026-May-05.xlsx"
)

# Generate timestamped output directory
timestamp   <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_root <- file.path(common_dir, paste0("bayesian_stats_histo_", timestamp))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Histology distance grid (used for factor levels and plot breaks)
distance_levels <- c(40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486)


######### Helper: normalize Region labels ####################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}


######### Helper: import one stain sheet, cleaned ############################
# Forces the response column to numeric (blank / text cells -> NA) and drops
# non-finite rows so the model fit and the raw means are built from identical,
# clean data on the same (z-score) scale.
import_stain <- function(input_file, sheet, response) {
  read_excel(input_file, sheet = sheet) %>%
    rename(!!response := stain_value, distance = Distance) %>%
    mutate(
      Region                = normalize_region(Region),
      Stage                 = as.integer(as.numeric(Stage)),
      Groups                = factor(Groups, levels = c("Control", "Experimental")),
      !!response            := suppressWarnings(as.numeric(.data[[response]]))
    ) %>%
    filter(is.finite(.data[[response]]))
}


######### Helper: fit one within-stage model #################################
# Within a fixed stain x region x source x stage, the only predictor is
# distance. Subject handling:
#   1 subject  -> response ~ distance_f                 (stan_glm)
#   2 subjects -> response ~ distance_f + subjectID     (stan_glm, fixed covar)
#   >=3        -> response ~ distance_f + (1|subjectID)  (stan_glmer, random int)
# The fitted curve is later averaged over whatever subjects are present.
fit_stage_model <- function(data, response, fit_ctrl) {
  n_subj <- length(unique(data$subjectID))
  
  if (n_subj >= 3) {
    f <- as.formula(paste0(response, " ~ distance_f + (1 | subjectID)"))
    stan_glmer(
      f, data = data, family = gaussian(),
      chains  = fit_ctrl$chains, iter   = fit_ctrl$iter,
      warmup  = fit_ctrl$warmup, cores  = fit_ctrl$cores,
      seed    = fit_ctrl$seed,   control = fit_ctrl$control
    )
  } else {
    f <- if (n_subj == 2) {
      as.formula(paste0(response, " ~ distance_f + subjectID"))
    } else {
      as.formula(paste0(response, " ~ distance_f"))
    }
    stan_glm(
      f, data = data, family = gaussian(),
      chains  = fit_ctrl$chains, iter   = fit_ctrl$iter,
      warmup  = fit_ctrl$warmup, cores  = fit_ctrl$cores,
      seed    = fit_ctrl$seed,   control = fit_ctrl$control
    )
  }
}


######### Helper: posterior fitted means by distance #########################
# Returns a draws x distance matrix of the model-implied mean at each distance,
# averaged over the subjects present in the fitting data. Uses only distances
# actually present in this stage's data so the prediction grid never references
# an unestimated factor level. If posterior_epred is unavailable in your
# rstanarm version, swap to:
#   ep <- posterior_linpred(model, newdata = newdata, transform = TRUE)
posterior_fitted_by_distance <- function(model, data) {
  d_present    <- sort(unique(data$distance))
  subj_present <- unique(data$subjectID)
  
  newdata <- expand.grid(
    distance  = d_present,
    subjectID = subj_present,
    stringsAsFactors = FALSE
  )
  newdata$distance_f <- factor(newdata$distance, levels = distance_levels)
  
  ep <- posterior_epred(model, newdata = newdata)   # draws x nrow(newdata)
  
  draws_by_d <- sapply(d_present, function(dd) {
    cols <- which(newdata$distance == dd)
    rowMeans(ep[, cols, drop = FALSE])
  })
  if (is.null(dim(draws_by_d))) draws_by_d <- matrix(draws_by_d, nrow = nrow(ep))
  colnames(draws_by_d) <- as.character(d_present)
  
  list(draws = draws_by_d, distances = d_present)
}


######### Method 1: model-based difference (stage k - stage 0) ###############
model_diff_table <- function(model0, data0, modelk, datak, prob = 0.95) {
  f0 <- posterior_fitted_by_distance(model0, data0)
  fk <- posterior_fitted_by_distance(modelk, datak)
  
  common <- sort(intersect(f0$distances, fk$distances))
  if (length(common) == 0) return(NULL)
  
  nd     <- min(nrow(f0$draws), nrow(fk$draws))   # align draw counts
  a_low  <- (1 - prob) / 2
  a_high <- 1 - a_low
  
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
# n counts only finite (non-blank) observations, so the propagated SE
#   se_diff = sqrt(s_k^2/n_k + s_0^2/n_0)
# reflects the true post-filter sample sizes. Raw components (s, n) are kept so
# a pooled SD sqrt(s_k^2 + s_0^2) can be recomputed if preferred.
raw_means_by_distance <- function(data, response) {
  data %>%
    group_by(distance) %>%
    summarise(
      mean = mean(.data[[response]], na.rm = TRUE),
      sd   = sd(.data[[response]],   na.rm = TRUE),
      n    = sum(is.finite(.data[[response]])),
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
    scale_x_continuous(breaks = distance_levels) +
    labs(
      x = "Distance", y = "Difference (Stage k - Stage 0)",
      title = title,
      subtitle = "blue = model fitted-curve difference (95% CI ribbon); red = raw-mean difference (\u00b1 SE)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}


######## BEGIN LOGIC #########################################################

fname      <- input_filenames[1]
input_file <- file.path(common_dir, fname)
message("ANALYZING DATA FOR: ", fname)

######## 1) IMPORT DATA ######################################################
# Groups encodes the measurement source:
#   Experimental = measurement around the enlarged perivascular space (EPVS)
#   Control      = measurement around the vessel
# Stage encodes CAA severity: 0 = control, 1 = mild, 2 = moderate, 3 = severe.
# stain_value is a z-score; blank / non-numeric cells are dropped on import.
lhe_data  <- import_stain(input_file, "lhe",  "lhe")
gfap_data <- import_stain(input_file, "gfap", "gfap")
cd68_data <- import_stain(input_file, "cd68", "cd68")

# Quick report of how many rows survived the blank/NA filter per sheet
message("Rows after filtering blanks -> lhe: ", nrow(lhe_data),
        " | gfap: ", nrow(gfap_data), " | cd68: ", nrow(cd68_data))

######## 2) DEFINE LOOP PARAMETERS ###########################################
stains         <- list(lhe = lhe_data, gfap = gfap_data, cd68 = cd68_data)
regions        <- c("front", "occip")
sources        <- c(epvs = "Experimental", vessel = "Control")   # name -> Groups level
nonzero_stages <- c(1L, 2L, 3L)
stage_labels   <- c("0" = "control", "1" = "mild", "2" = "moderate", "3" = "severe")

fit_ctrl <- list(
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,          # must be < iter; yields 2000 post-warmup draws/chain
  cores   = 4,
  seed    = 42,
  control = list(
    adapt_delta   = 0.95,
    max_treedepth = 15
  )
)
base_name <- tools::file_path_sans_ext(fname)

######## 3) MAIN LOOP: STAIN x REGION x SOURCE x (STAGE k - STAGE 0) #########
for (stain in names(stains)) {
  sdata <- stains[[stain]]
  
  for (rgn in regions) {
    for (src in names(sources)) {
      grp <- sources[[src]]
      
      message("\n===== ", stain, " | ", rgn, " | ", src, " (Groups = ", grp, ") =====")
      
      base_sub <- sdata %>%
        filter(Region == rgn, Groups == grp) %>%
        mutate(distance_f = factor(distance, levels = distance_levels)) %>%
        droplevels()
      
      data0 <- base_sub %>% filter(Stage == 0)
      if (nrow(data0) == 0) {
        message("  Skipping - no stage 0 (control) data for this stain/region/source")
        next
      }
      
      # Fit the stage 0 model once and reuse for every nonzero-stage comparison
      message("  Fitting stage 0 (control) model...")
      model0 <- tryCatch(
        fit_stage_model(data0, stain, fit_ctrl),
        error = function(e) { message("  stage 0 fit FAILED: ", e$message); NULL }
      )
      if (!inherits(model0, "stanreg")) {
        message("  Skipping - stage 0 model did not fit for this stain/region/source")
        next
      }
      
      out_dir <- file.path(output_root, "stage_vs_control", rgn, src)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      model_accum <- list()
      raw_accum   <- list()
      
      for (k in nonzero_stages) {
        datak <- base_sub %>% filter(Stage == k)
        if (nrow(datak) == 0) {
          message("  Skipping stage ", k, " (", stage_labels[as.character(k)],
                  ") - no data for this stain/region/source")
          next
        }
        
        message("  Fitting stage ", k, " (", stage_labels[as.character(k)], ") model...")
        modelk <- tryCatch(
          fit_stage_model(datak, stain, fit_ctrl),
          error = function(e) { message("  stage ", k, " fit FAILED: ", e$message); NULL }
        )
        if (!inherits(modelk, "stanreg")) {
          message("  Skipping stage ", k, " - model did not fit")
          next
        }
        
        ## Method 1: model fitted-curve difference
        mdiff <- model_diff_table(model0, data0, modelk, datak)
        ## Method 2: raw-mean difference (SE recomputed on post-filter n)
        rdiff <- raw_diff_table(data0, datak, stain)
        
        if (is.null(mdiff) || nrow(rdiff) == 0) {
          message("  No common distances between stage ", k, " and stage 0 - skipping")
          next
        }
        
        ## Save both summary tables
        model_fout <- file.path(
          out_dir,
          paste0(base_name, "__", stain, "_", rgn, "_", src,
                 "_stage", k, "_vs_stage0_model_diff.csv")
        )
        raw_fout <- file.path(
          out_dir,
          paste0(base_name, "__", stain, "_", rgn, "_", src,
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
          title = paste0(stain, " (", rgn, ", ", src, "): stage - control")
        )
        print(p)
      }
      
      message("  Finished: ", out_dir)
    }
  }
}