# load necessary libraries
library(rstanarm)
library(tidyverse)   # dplyr, tidyr, stringr, ggplot2, purrr, readr


######## Filepaths ###########################################################
common_dir <- "/projectnb/npbssmic/ns/CAA"
# Output of measure_parenchyma_shells.m + add_percent_change_to_csv.m
input_filenames <- c(
  "parenchyma_shells_2026-10-07_pct.csv"
)


######## Settings ############################################################
# Label each 40 µm shell by its outer edge: 0-40 -> 40, ..., 480-520 -> 520.
# Use "d_inner_um" instead to label by inner edge (0 ... 480).
distance_col <- "d_outer_um"

stage_levels   <- c("control", "moderate", "severe")
nonzero_stages <- c("moderate", "severe")   # each compared against control
stage_colors   <- c(control = "grey35", moderate = "darkorange", severe = "firebrick")

# MATLAB class -> source name used in output folders
sources <- c(epvs = "EPVS", vessel = "vessel")

# Two separate model sets: raw values and percentage change from WM median.
# CSV column for each property x scale:
response_cols <- list(
  scattering = c(raw = "mus", pct = "mus_pct"),
  retardance = c(raw = "ret", pct = "ret_pct")
)
scales      <- c("raw", "pct")
scale_ylab  <- list(
  raw = c(scattering = "Scattering coefficient (raw)",
          retardance = "Retardance (raw)"),
  pct = c(scattering = "Scattering coefficient (% change from WM median)",
          retardance = "Retardance (% change from WM median)")
)

# Raw-value bounds (mus_ok / ret_ok from the CSV). Applied to BOTH scales so
# the raw and percentage models are fit to exactly the same rows.
apply_bounds <- TRUE
ok_cols      <- c(scattering = "mus_ok", retardance = "ret_ok")

# stan_glmer with one intercept per object is slower than stan_glm;
# drop iter to 2000 if runtime is a problem and diagnostics stay clean.
fit_ctrl <- list(chains = 4, iter = 4000, cores = 4, seed = 42, adapt_delta = 0.95)

output_filepath <- "/projectnb/npbssmic/ns/CAA/beta_stats"
save_models     <- TRUE
# Resume: if a model's .rds already exists, load it instead of refitting.
# Delete the saved .rds files if the data or the model formula change.
reuse_saved_models <- TRUE

# Figure export
save_plots   <- TRUE
plot_formats <- c("png", "jpeg")   # any ggsave device: "png", "jpeg", "pdf", "tiff"
plot_dpi     <- 300
plot_width   <- 7                  # inches
plot_height  <- 5

# Full posterior draws (rows = draws, columns = distances) for replotting
save_draws <- TRUE


######### Helper: normalize Region labels ####################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}


######### Helper: fit one within-stage model #################################
# Within a fixed optical property x region x source x stage, distance is the
# predictor of interest. Each object (vessel or EPVS) contributes one value per
# shell, so (1 | object_id) accounts for repeated shells from the same object.
# subjectID enters as a fixed covariate when the stage has >= 2 subjects
# (control, severe); moderate has one subject, so the term is dropped.
fit_stage_model <- function(data, response, fit_ctrl) {
  n_subj <- n_distinct(data$subjectID)
  rhs <- if (n_subj >= 2) "distance_f + subjectID" else "distance_f"
  f <- as.formula(paste0(response, " ~ ", rhs, " + (1 | object_id)"))
  stan_glmer(
    f, data = data, family = gaussian(),
    chains = fit_ctrl$chains, iter = fit_ctrl$iter,
    cores  = fit_ctrl$cores,  seed = fit_ctrl$seed,
    adapt_delta = fit_ctrl$adapt_delta
  )
}


######### Helper: fit (or reload) a model and report convergence #############
report_diagnostics <- function(m, label) {
  sm   <- summary(m)
  rhat <- max(sm[, "Rhat"],  na.rm = TRUE)
  ness <- min(sm[, "n_eff"], na.rm = TRUE)
  worst <- rownames(sm)[which.min(sm[, "n_eff"])]
  flag <- if (rhat > 1.01 || ness < 400) "  <-- CHECK" else ""
  message(sprintf("  %s diagnostics: max Rhat = %.3f, min n_eff = %.0f (%s)%s",
                  label, rhat, ness, worst, flag))
}

fit_or_load <- function(data, response, fit_ctrl, rds_path, label, legacy_path = NULL) {
  if (reuse_saved_models && file.exists(rds_path)) {
    message("  Loading saved ", label, " model: ", basename(rds_path))
    m <- readRDS(rds_path)
  } else if (reuse_saved_models && !is.null(legacy_path) && file.exists(legacy_path)) {
    # Model saved by an earlier run in the old flat folder layout: reuse it
    # and copy it into the new stage subfolder
    message("  Loading saved ", label, " model (old location): ", basename(legacy_path))
    m <- readRDS(legacy_path)
    if (save_models) saveRDS(m, rds_path)
  } else {
    m <- fit_stage_model(data, response, fit_ctrl)
    if (save_models) saveRDS(m, rds_path)
  }
  report_diagnostics(m, label)
  m
}


######### Helper: posterior fitted means by distance #########################
# Returns a draws x distance matrix of the model-implied mean at each distance
# for a typical object (object random effect = 0, via re.form = NA), averaged
# over the subjects present in the fitting data. Only distances present in
# this stage's data are predicted.
posterior_fitted_by_distance <- function(model, data) {
  d_present    <- sort(unique(data$distance))
  subj_present <- unique(data$subjectID)
  
  newdata <- expand.grid(
    distance  = d_present,
    subjectID = subj_present,
    stringsAsFactors = FALSE
  )
  newdata$distance_f <- factor(newdata$distance, levels = levels(data$distance_f))
  
  ep <- posterior_epred(model, newdata = newdata, re.form = NA)  # draws x nrow(newdata)
  
  draws_by_d <- sapply(d_present, function(dd) {
    cols <- which(newdata$distance == dd)
    rowMeans(ep[, cols, drop = FALSE])
  })
  if (is.null(dim(draws_by_d))) draws_by_d <- matrix(draws_by_d, nrow = nrow(ep))
  colnames(draws_by_d) <- as.character(d_present)
  
  list(draws = draws_by_d, distances = d_present)
}


######### Per-stage posterior summary (absolute fitted curve) ################
# NOTE ON INTERPRETATION: subjectID is a fixed covariate (or absent for
# moderate), so the interval describes the measured subjects only. The object
# random intercept accounts for repeated shells within an object; it does not
# propagate between-subject variance.
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


######### Method 1: model-based difference (stage k - control) ###############
model_diff_table <- function(model0, data0, modelk, datak, prob = 0.95) {
  f0 <- posterior_fitted_by_distance(model0, data0)
  fk <- posterior_fitted_by_distance(modelk, datak)
  
  common <- sort(intersect(f0$distances, fk$distances))
  if (length(common) == 0) return(NULL)
  
  nd     <- min(nrow(f0$draws), nrow(fk$draws))
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


######### Method 2: raw-mean difference (stage k - control) ##################
# Each row is one object at one distance, so n = number of objects contributing
# to that distance.
raw_means_by_distance <- function(data, response) {
  data %>%
    group_by(distance) %>%
    summarise(
      mean = mean(.data[[response]], na.rm = TRUE),
      sd   = sd(.data[[response]],   na.rm = TRUE),
      n    = sum(!is.na(.data[[response]])),
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
plot_stage_curves <- function(fit_df, raw_df, title, ylab, breaks) {
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
    scale_x_continuous(breaks = breaks) +
    scale_color_manual(values = stage_colors, name = NULL,
                       aesthetics = c("color", "fill")) +
    labs(
      x = "Distance from boundary (\u00b5m, shell outer edge)", y = ylab,
      title = title,
      subtitle = "lines/ribbons = posterior fitted mean (95% CI); points = raw means across objects (\u00b1 SE)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          legend.position = "top")
}


######### Plot: overlay of both difference methods ###########################
plot_stage_diffs <- function(model_df, raw_df, title, breaks,
                             ylab = "Difference (stage - control)") {
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
    scale_x_continuous(breaks = breaks) +
    labs(
      x = "Distance from boundary (\u00b5m, shell outer edge)",
      y = ylab,
      title = title,
      subtitle = "blue = model fitted-curve difference (95% CI ribbon); red = raw-mean difference (\u00b1 SE)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8))
}


######### Helpers: export figures and posterior draws ########################
save_plot <- function(p, stem_path, width = plot_width, height = plot_height) {
  if (!save_plots) return(invisible(NULL))
  for (fmt in plot_formats) {
    ggsave(paste0(stem_path, ".", fmt), p, width = width, height = height,
           dpi = plot_dpi, bg = "white")   # white bg: JPEG has no transparency
  }
}

# Wide CSV: column 'draw', then one column per distance (d40, d80, ...).
# MATLAB: T = readtable(file); M = T{:, 2:end};
write_draws <- function(draws, path) {
  if (!save_draws) return(invisible(NULL))
  out <- as.data.frame(draws)
  names(out) <- paste0("d", colnames(draws))
  out <- cbind(draw = seq_len(nrow(out)), out)
  write.csv(out, path, row.names = FALSE)
}

diff_draws <- function(model0, data0, modelk, datak) {
  f0 <- posterior_fitted_by_distance(model0, data0)
  fk <- posterior_fitted_by_distance(modelk, datak)
  common <- as.character(sort(intersect(f0$distances, fk$distances)))
  nd <- min(nrow(f0$draws), nrow(fk$draws))
  fk$draws[seq_len(nd), common, drop = FALSE] - f0$draws[seq_len(nd), common, drop = FALSE]
}


######## BEGIN LOGIC #########################################################

fname      <- input_filenames[1]
input_file <- file.path(common_dir, fname)
message("ANALYZING DATA FOR: ", fname)

######## 1) IMPORT DATA ######################################################
# One row per object (vessel or EPVS) x 40 µm shell.
#   class    = "EPVS" or "vessel" (measurement source)
#   severity = control / moderate / severe
raw <- read_csv(input_file, show_col_types = FALSE)

# MATLAB may write logicals as 1/0 or true/false
as_flag <- function(x) as.character(x) %in% c("TRUE", "true", "1")

needed <- c("subject", "severity", "region", "class", "object_id", distance_col,
            "enough_voxels", unlist(response_cols), ok_cols)
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) stop("CSV is missing columns: ", paste(missing_cols, collapse = ", "))

keep_row <- as_flag(raw$enough_voxels)
message("Rows with enough voxels: ", sum(keep_row), " of ", nrow(raw))

shells <- raw[keep_row, ] %>%
  transmute(
    subjectID  = as.character(subject),
    Stage      = factor(severity, levels = stage_levels),
    Region     = normalize_region(region),
    Groups     = as.character(class),
    object_id  = as.character(object_id),
    distance   = .data[[distance_col]],
    n_vox, obj_diam_um,
    mus_ok     = as_flag(mus_ok),
    ret_ok     = as_flag(ret_ok),
    scattering_raw = mus,
    scattering_pct = mus_pct,
    retardance_raw = ret,
    retardance_pct = ret_pct
  )

# Apply raw-value bounds to both scales of each property
if (apply_bounds) {
  message("Rows failing bounds: scattering ", sum(!shells$mus_ok),
          ", retardance ", sum(!shells$ret_ok))
  shells <- shells %>%
    mutate(
      across(c(scattering_raw, scattering_pct), ~ if_else(mus_ok, .x, NA_real_)),
      across(c(retardance_raw, retardance_pct), ~ if_else(ret_ok, .x, NA_real_))
    )
}

## Guards
bad_stage <- setdiff(unique(raw$severity), stage_levels)
if (length(bad_stage) > 0) {
  warning("Unexpected severity labels (rows ignored): ", paste(bad_stage, collapse = ", "))
  shells <- shells %>% filter(!is.na(Stage))
}
dist_levels <- sort(unique(shells$distance))
message("Distances present: ", paste(dist_levels, collapse = ", "))
print(shells %>% count(Stage, Region, Groups, subjectID, name = "rows") %>%
        left_join(shells %>% group_by(Stage, Region, Groups, subjectID) %>%
                    summarise(objects = n_distinct(object_id), .groups = "drop"),
                  by = c("Stage", "Region", "Groups", "subjectID")), n = Inf)

######## 2) DEFINE LOOP PARAMETERS ###########################################
properties <- names(response_cols)
regions    <- c("front", "occip")
base_name  <- tools::file_path_sans_ext(fname)

prep_stage <- function(base_sub, stage) {
  base_sub %>%
    filter(Stage == stage) %>%
    mutate(distance_f = factor(distance, levels = dist_levels)) %>%
    droplevels()            # no coefficients for distances absent in this stage
}

######## 3) MAIN LOOP: PROPERTY x REGION x SOURCE x (STAGE - CONTROL) ########
for (scale in scales) {
  for (prop in properties) {
    resp <- paste0(prop, "_", scale)          # response column in `shells`
    ylab <- scale_ylab[[scale]][[prop]]
    for (rgn in regions) {
      for (src in names(sources)) {
        grp <- sources[[src]]
        
        message("\n===== ", scale, " | ", prop, " | ", rgn, " | ", src, " =====")
        
        base_sub <- shells %>%
          filter(Region == rgn, Groups == grp, !is.na(.data[[resp]]))
        
        data0 <- prep_stage(base_sub, "control")
        if (nrow(data0) == 0) {
          message("  Skipping - no control data for this property/region/source")
          next
        }
        
        # Layout: <block_dir>/control/, /severe/, /moderate/,
        #         /severe_vs_control/, /moderate_vs_control/, /all_stages/
        block_dir <- file.path(output_filepath, base_name, scale, "stage_vs_control", rgn, src)
        sub_dir <- function(name) {
          d <- file.path(block_dir, name)
          dir.create(d, recursive = TRUE, showWarnings = FALSE)
          d
        }
        ctrl_dir <- sub_dir("control")
        all_dir  <- sub_dir("all_stages")
        stem <- paste0(base_name, "__", scale, "_", prop, "_", rgn, "_", src)
        
        message("  Fitting control model: ", n_distinct(data0$subjectID), " subject(s), ",
                n_distinct(data0$object_id), " objects, ", nrow(data0), " rows")
        model0 <- fit_or_load(data0, resp, fit_ctrl,
                              file.path(ctrl_dir, paste0(stem, "_control_model.rds")), "control",
                              legacy_path = file.path(block_dir, paste0(stem, "_control_model.rds")))
        
        write_draws(posterior_fitted_by_distance(model0, data0)$draws,
                    file.path(ctrl_dir, paste0(stem, "_control_draws.csv")))
        
        fit0 <- posterior_stage_table(model0, data0) %>%
          mutate(stage_label = "control", .before = 1)
        raw0 <- raw_means_by_distance(data0, resp) %>%
          mutate(stage_label = "control", .before = 1)
        write.csv(fit0, file.path(ctrl_dir, paste0(stem, "_control_fitted.csv")), row.names = FALSE)
        write.csv(raw0, file.path(ctrl_dir, paste0(stem, "_control_raw.csv")),    row.names = FALSE)
        
        fit_accum   <- list(control = fit0)
        raw_accum_s <- list(control = raw0)
        model_accum <- list()
        raw_accum   <- list()
        
        for (k in nonzero_stages) {
          datak <- prep_stage(base_sub, k)
          if (nrow(datak) == 0) {
            message("  Skipping ", k, " - no data for this property/region/source")
            next
          }
          
          stage_dir <- sub_dir(k)
          diff_dir  <- sub_dir(paste0(k, "_vs_control"))
          
          message("  Fitting ", k, " model: ", n_distinct(datak$subjectID), " subject(s), ",
                  n_distinct(datak$object_id), " objects, ", nrow(datak), " rows")
          modelk <- fit_or_load(datak, resp, fit_ctrl,
                                file.path(stage_dir, paste0(stem, "_", k, "_model.rds")), k,
                                legacy_path = file.path(block_dir, paste0(stem, "_", k, "_model.rds")))
          
          write_draws(posterior_fitted_by_distance(modelk, datak)$draws,
                      file.path(stage_dir, paste0(stem, "_", k, "_draws.csv")))
          
          fitk <- posterior_stage_table(modelk, datak) %>%
            mutate(stage_label = k, .before = 1)
          rawk <- raw_means_by_distance(datak, resp) %>%
            mutate(stage_label = k, .before = 1)
          write.csv(fitk, file.path(stage_dir, paste0(stem, "_", k, "_fitted.csv")), row.names = FALSE)
          write.csv(rawk, file.path(stage_dir, paste0(stem, "_", k, "_raw.csv")),    row.names = FALSE)
          
          fit_accum[[k]]   <- fitk
          raw_accum_s[[k]] <- rawk
          
          mdiff <- model_diff_table(model0, data0, modelk, datak)
          rdiff <- raw_diff_table(data0, datak, resp)
          if (is.null(mdiff) || nrow(rdiff) == 0) {
            message("  No common distances between ", k, " and control - skipping contrast")
            next
          }
          write.csv(mdiff, file.path(diff_dir, paste0(stem, "_", k, "_vs_control_model_diff.csv")),
                    row.names = FALSE)
          write.csv(rdiff, file.path(diff_dir, paste0(stem, "_", k, "_vs_control_raw_diff.csv")),
                    row.names = FALSE)
          write_draws(diff_draws(model0, data0, modelk, datak),
                      file.path(diff_dir, paste0(stem, "_", k, "_vs_control_diff_draws.csv")))
          
          sl <- paste0(k, " - control")
          model_accum[[k]] <- mdiff %>% mutate(stage_label = sl)
          raw_accum[[k]]   <- rdiff %>% mutate(stage_label = sl)
        }
        
        fit_all       <- bind_rows(fit_accum)     %>% mutate(stage_label = factor(stage_label, stage_levels))
        raw_all_stage <- bind_rows(raw_accum_s)   %>% mutate(stage_label = factor(stage_label, stage_levels))
        write.csv(fit_all,       file.path(all_dir, paste0(stem, "_all_stages_fitted.csv")), row.names = FALSE)
        write.csv(raw_all_stage, file.path(all_dir, paste0(stem, "_all_stages_raw.csv")),    row.names = FALSE)
        
        p_curves <- plot_stage_curves(
          fit_all, raw_all_stage,
          title  = paste0(prop, " [", scale, "] (", rgn, ", ", src, "): fitted curve by stage"),
          ylab   = ylab, breaks = dist_levels
        )
        print(p_curves)
        save_plot(p_curves, file.path(all_dir, paste0(stem, "_curves")))
        
        if (length(model_accum) > 0) {
          model_all <- bind_rows(model_accum)
          raw_all   <- bind_rows(raw_accum)
          write.csv(model_all, file.path(all_dir, paste0(stem, "_all_stages_model_diff.csv")),
                    row.names = FALSE)
          write.csv(raw_all,   file.path(all_dir, paste0(stem, "_all_stages_raw_diff.csv")),
                    row.names = FALSE)
          p_diff <- plot_stage_diffs(
            model_all, raw_all,
            title  = paste0(prop, " [", scale, "] (", rgn, ", ", src, "): stage - control"),
            ylab   = paste0("Difference (stage - control), ", ylab),
            breaks = dist_levels
          )
          print(p_diff)
          save_plot(p_diff, file.path(all_dir, paste0(stem, "_diffs")),
                    width = plot_width * max(1, length(model_accum)) * 0.8)
        }
        
        message("  Finished: ", block_dir)
      }
    }
  }
}  # scale