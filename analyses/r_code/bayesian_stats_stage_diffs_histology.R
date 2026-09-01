# =============================================================================
# Bayesian per-stage fitted curves + severity contrasts for HISTOPATHOLOGY
# (CD68, GFAP, LHE), modeled over annular-ring distance bins around EPVS/vessel.
#
# This script now mirrors the *function* of the optical-property script
# (retardance / scattering). In addition to the stage-k - stage-0 contrasts it
# already produced, it emits, for every stain x region x source x stage:
#   - the absolute posterior fitted curve (mean, SD, 95% CI) by distance
#   - the absolute raw means (+/- SE) by distance
#   - combined "all stages" fitted/raw tables
#   - an absolute-curve plot (one line/ribbon per stage)
#   - optionally, the saved stanreg objects (save_models)
# alongside the existing model/raw difference tables and the difference plot.
# =============================================================================

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

save_models <- TRUE   # keep stanreg objects so summaries can be recomputed
# without refitting


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


######### Per-stage posterior summary (absolute fitted curve) ################
# Summarizes a single stage's fitted curve rather than a contrast: posterior
# mean, SD and central interval of the model-implied mean at each distance,
# averaged over the subjects present in that stage.
#
# NOTE ON INTERPRETATION: for the 1- and 2-subject stages subjectID is a fixed
# covariate, so this interval describes the mean of the subjects actually
# measured. It does not propagate between-subject variance and should not be
# read as generalizing to new subjects at this stage. For >= 3 subjects the
# random-intercept fit is likewise averaged over the measured subjects' fitted
# values (subject random effects included via posterior_epred).
posterior_stage_table <- function(model, data, prob = 0.95) {
  f <- posterior_fitted_by_distance(model, data)
  a_low  <- (1 - prob) / 2
  a_high <- 1 - a_low
  
  purrr::map_dfr(seq_along(f$distances), function(i) {
    v <- f$draws[, i]
    tibble(
      distance = f$distances[i],
      fit_mean = mean(v),
      fit_sd   = sd(v),
      lower    = unname(quantile(v, a_low)),
      upper    = unname(quantile(v, a_high)),
      n_draws  = length(v)
    )
  }) %>% arrange(distance)
}


######### Method 1: model-based difference (stage k - stage 0) ###############
# Differences the two independent stage models' posterior fitted curves draw by
# draw, then summarizes the difference distribution at each common distance.
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
# a pooled SD sqrt(s_k^2 + s_0^2) can be recomputed if preferred. An `se` column
# (sd/sqrt(n)) is also returned for the absolute per-stage error bars; the
# difference table selects it away to avoid a join-side name collision.
raw_means_by_distance <- function(data, response) {
  data %>%
    group_by(distance) %>%
    summarise(
      mean = mean(.data[[response]], na.rm = TRUE),
      sd   = sd(.data[[response]],   na.rm = TRUE),
      n    = sum(is.finite(.data[[response]])),
      .groups = "drop"
    ) %>%
    mutate(se = sd / sqrt(n))
}

raw_diff_table <- function(data0, datak, response) {
  r0 <- raw_means_by_distance(data0, response) %>%
    select(distance, mean_0 = mean, sd_0 = sd, n_0 = n)
  rk <- raw_means_by_distance(datak, response) %>%
    select(distance, mean_k = mean, sd_k = sd, n_k = n)
  
  inner_join(r0, rk, by = "distance") %>%
    mutate(
      diff    = mean_k - mean_0,
      se_diff = sqrt(sd_k^2 / n_k + sd_0^2 / n_0)
    ) %>%
    arrange(distance) %>%
    select(distance, mean_0, mean_k, diff, sd_0, sd_k, n_0, n_k, se_diff)
}


######### Plot: absolute fitted curves, one line per stage ###################
# Two stages in the corrected design; scale_manual maps each label to a colour.
stage_palette <- c(
  control = "grey35",
  severe  = "firebrick"
)

plot_stage_curves <- function(fit_df, raw_df, title, ylab,
                              palette = stage_palette) {
  ggplot() +
    geom_ribbon(data = fit_df,
                aes(x = distance, ymin = lower, ymax = upper, fill = stage_label),
                alpha = 0.18) +
    geom_line(data = fit_df,
              aes(x = distance, y = fit_mean, color = stage_label), linewidth = 0.8) +
    geom_point(data = raw_df,
               aes(x = distance, y = mean, color = stage_label), size = 1.2) +
    geom_errorbar(data = raw_df,
                  aes(x = distance, ymin = mean - se, ymax = mean + se,
                      color = stage_label),
                  width = 8, alpha = 0.6) +
    scale_x_continuous(breaks = distance_levels) +
    scale_color_manual(values = palette, name = NULL,
                       aesthetics = c("color", "fill")) +
    labs(
      x = "Distance", y = ylab,
      title = title,
      subtitle = "lines/ribbons = posterior fitted mean (95% CI); points = raw means (\u00b1 SE)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}


######### Plot: overlay of both difference methods ###########################
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
# Stage encodes CAA severity. The corrected design retains only:
#   0 = control (n = 2 subjects), 3 = severe (n = 2 subjects).
# stain_value is a z-score; blank / non-numeric cells are dropped on import.
lhe_data  <- import_stain(input_file, "lhe",  "lhe")
gfap_data <- import_stain(input_file, "gfap", "gfap")
cd68_data <- import_stain(input_file, "cd68", "cd68")

# Quick report of how many rows survived the blank/NA filter per sheet
message("Rows after filtering blanks -> lhe: ", nrow(lhe_data),
        " | gfap: ", nrow(gfap_data), " | cd68: ", nrow(cd68_data))

## Guard: confirm the workbook really contains only stages 0 and 3
observed_stages <- sort(unique(c(lhe_data$Stage, gfap_data$Stage, cd68_data$Stage)))
message("Stages present in input: ", paste(observed_stages, collapse = ", "))
if (!all(observed_stages %in% c(0L, 3L))) {
  warning("Unexpected stage codes present: ",
          paste(setdiff(observed_stages, c(0L, 3L)), collapse = ", "),
          " -- these rows will be ignored.")
}

######## 2) DEFINE LOOP PARAMETERS ###########################################
stains         <- list(lhe = lhe_data, gfap = gfap_data, cd68 = cd68_data)
regions        <- c("front", "occip")
sources        <- c(epvs = "Experimental", vessel = "Control")   # name -> Groups level
# Corrected 2x2 design: only control (0) and severe (3), so the sole contrast
# is severe - control.
nonzero_stages <- c(3L)
stage_labels   <- c("0" = "control", "3" = "severe")

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
      
      out_dir <- file.path(output_root, "stage_vs_control", rgn, src)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      stem <- paste0(base_name, "__", stain, "_", rgn, "_", src)
      
      # Fit the stage 0 model once and reuse for every nonzero-stage comparison
      message("  Fitting stage 0 (control) model on ",
              length(unique(data0$subjectID)), " subject(s)...")
      model0 <- tryCatch(
        fit_stage_model(data0, stain, fit_ctrl),
        error = function(e) { message("  stage 0 fit FAILED: ", e$message); NULL }
      )
      if (!inherits(model0, "stanreg")) {
        message("  Skipping - stage 0 model did not fit for this stain/region/source")
        next
      }
      if (save_models) {
        saveRDS(model0, file.path(out_dir, paste0(stem, "_stage0_model.rds")))
      }
      
      ## --- Per-stage summaries for stage 0 -------------------------------
      fit0 <- posterior_stage_table(model0, data0) %>%
        mutate(stage = 0L, stage_label = stage_labels["0"], .before = 1)
      raw0 <- raw_means_by_distance(data0, stain) %>%
        mutate(stage = 0L, stage_label = stage_labels["0"], .before = 1)
      
      write.csv(fit0, file.path(out_dir, paste0(stem, "_stage0_fitted.csv")),
                row.names = FALSE)
      write.csv(raw0, file.path(out_dir, paste0(stem, "_stage0_raw.csv")),
                row.names = FALSE)
      
      message("  --- Stage 0 (control) posterior fitted curve ---")
      print(fit0, n = Inf)
      
      fit_accum   <- list("0" = fit0)
      raw_accum_s <- list("0" = raw0)
      model_accum <- list()
      raw_accum   <- list()
      
      for (k in nonzero_stages) {
        datak <- base_sub %>% filter(Stage == k)
        if (nrow(datak) == 0) {
          message("  Skipping stage ", k, " (", stage_labels[as.character(k)],
                  ") - no data for this stain/region/source")
          next
        }
        
        message("  Fitting stage ", k, " (", stage_labels[as.character(k)],
                ") model on ", length(unique(datak$subjectID)), " subject(s)...")
        modelk <- tryCatch(
          fit_stage_model(datak, stain, fit_ctrl),
          error = function(e) { message("  stage ", k, " fit FAILED: ", e$message); NULL }
        )
        if (!inherits(modelk, "stanreg")) {
          message("  Skipping stage ", k, " - model did not fit")
          next
        }
        if (save_models) {
          saveRDS(modelk, file.path(out_dir, paste0(stem, "_stage", k, "_model.rds")))
        }
        
        ## --- Per-stage summaries for stage k -----------------------------
        fitk <- posterior_stage_table(modelk, datak) %>%
          mutate(stage = k, stage_label = stage_labels[as.character(k)], .before = 1)
        rawk <- raw_means_by_distance(datak, stain) %>%
          mutate(stage = k, stage_label = stage_labels[as.character(k)], .before = 1)
        
        write.csv(fitk, file.path(out_dir, paste0(stem, "_stage", k, "_fitted.csv")),
                  row.names = FALSE)
        write.csv(rawk, file.path(out_dir, paste0(stem, "_stage", k, "_raw.csv")),
                  row.names = FALSE)
        
        message("  --- Stage ", k, " (", stage_labels[as.character(k)],
                ") posterior fitted curve ---")
        print(fitk, n = Inf)
        
        fit_accum[[as.character(k)]]   <- fitk
        raw_accum_s[[as.character(k)]] <- rawk
        
        ## --- Contrasts ---------------------------------------------------
        ## Method 1: model fitted-curve difference
        mdiff <- model_diff_table(model0, data0, modelk, datak)
        ## Method 2: raw-mean difference (SE recomputed on post-filter n)
        rdiff <- raw_diff_table(data0, datak, stain)
        
        if (is.null(mdiff) || nrow(rdiff) == 0) {
          message("  No common distances between stage ", k, " and stage 0 - skipping")
          next
        }
        
        write.csv(mdiff,
                  file.path(out_dir, paste0(stem, "_stage", k, "_vs_stage0_model_diff.csv")),
                  row.names = FALSE)
        write.csv(rdiff,
                  file.path(out_dir, paste0(stem, "_stage", k, "_vs_stage0_raw_diff.csv")),
                  row.names = FALSE)
        
        message("  --- Method 1 (model) : stage ", k, " - stage 0 ---")
        print(mdiff, n = Inf)
        message("  --- Method 2 (raw)   : stage ", k, " - stage 0 ---")
        print(rdiff, n = Inf)
        
        sl <- paste0("Stage ", k, " (", stage_labels[as.character(k)], ") - control")
        model_accum[[as.character(k)]] <- mdiff %>% mutate(stage = k, stage_label = sl)
        raw_accum[[as.character(k)]]   <- rdiff %>% mutate(stage = k, stage_label = sl)
      }
      
      ## --- Combined per-stage tables (all stages stacked) ----------------
      fit_all       <- bind_rows(fit_accum)
      raw_all_stage <- bind_rows(raw_accum_s)
      write.csv(fit_all,
                file.path(out_dir, paste0(stem, "_all_stages_fitted.csv")),
                row.names = FALSE)
      write.csv(raw_all_stage,
                file.path(out_dir, paste0(stem, "_all_stages_raw.csv")),
                row.names = FALSE)
      
      ## --- Absolute-curve plot -------------------------------------------
      p_curves <- plot_stage_curves(
        fit_all, raw_all_stage,
        title = paste0(stain, " (", rgn, ", ", src, "): fitted curve by stage"),
        ylab  = stain
      )
      print(p_curves)
      
      ## --- Difference plot -----------------------------------------------
      if (length(model_accum) > 0) {
        model_all <- bind_rows(model_accum)
        raw_all   <- bind_rows(raw_accum)
        p_diff <- plot_stage_diffs(
          model_all, raw_all,
          title = paste0(stain, " (", rgn, ", ", src, "): stage - control")
        )
        print(p_diff)
      }
      
      message("  Finished: ", out_dir)
    }
  }
}