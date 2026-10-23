###############################################################################
# Bayesian analysis of parenchymal optical properties around EPVS and vessels
# (brms version)
#
# For each scale x property x region:
#   A) GROUP, SEPARATE : EPVS and vessel fit in separate models, one per stage
#   B) GROUP, JOINT    : EPVS + vessel in one model per stage
#                        -> EPVS - vessel at each distance, and
#                           (EPVS - vessel)_stage - (EPVS - vessel)_control
#   C) SUBJECT         : one joint model per subject
#                        -> per-subject class curves and EPVS - vessel
#
# Input: output of measure_parenchyma_shells.m + add_percent_change_to_csv.m
###############################################################################

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(tidyverse)   # dplyr, tidyr, stringr, ggplot2, purrr, readr
})


######## Filepaths ###########################################################
common_dir      <- "/projectnb/npbssmic/ns/CAA"
input_filename  <- "parenchyma_shells_2026-10-07_pct.csv"
output_filepath <- "/projectnb/npbssmic/ns/CAA/beta_stats"


######## Analysis settings ###################################################
# Label each 40 µm shell by its outer edge (0-40 -> 40, ..., 480-520 -> 520)
distance_col <- "d_outer_um"

stage_levels   <- c("control", "moderate", "severe")
nonzero_stages <- c("moderate", "severe")          # each compared with control
stage_colors   <- c(control = "grey35", moderate = "darkorange", severe = "firebrick")

class_levels <- c("vessel", "EPVS")               # vessel = reference level
class_colors <- c(vessel = "grey35", EPVS = "steelblue")
sources      <- c(epvs = "EPVS", vessel = "vessel")

# Define the dataset to use
response_cols <- list(
  scattering = c(raw = "mus", pct = "mus_pct"),
  retardance = c(raw = "ret", pct = "ret_pct")
)
scales     <- c("pct") # c("raw", "pct") for both
regions    <- c("front", "occip")
scale_ylab <- list(
  raw = c(scattering = "Scattering coefficient (raw)",
          retardance = "Retardance (raw)"),
  pct = c(scattering = "Scattering coefficient (% change from WM median)",
          retardance = "Retardance (% change from WM median)")
)

# Raw-value bounds (mus_ok / ret_ok), applied to both scales
apply_bounds <- TRUE

# Which analyses to run
run_separate <- TRUE    # A) EPVS and vessel in separate models
run_joint    <- TRUE    # B) EPVS + vessel in one model per stage
run_subject  <- TRUE    # C) one joint model per subject

# Vessel subsampling: at most this many vessel objects per subject x region
# (all EPVS are always kept). NULL = use every vessel.
# vessel_subsample_n <- 500
vessel_subsample_n <- NULL
subsample_seed     <- 2026

# Diameter sensitivity analysis:
#   "none"      main analysis
#   "covariate" add standardized log(object diameter) as a fixed effect;
#               curves are evaluated at the mean diameter
#   "matched"   keep only vessels whose diameter lies within the EPVS
#               diameter range (matched_quantiles) of the same region
diameter_mode     <- "none"
matched_quantiles <- c(0.05, 0.95)

# Distances supported by fewer objects than this are flagged (hollow points)
min_objects_flag <- 10


######## brms settings #######################################################
model_family    <- "gaussian"   # or "student" for heavy-tailed residuals
heteroscedastic <- TRUE         # sigma ~ distance (+ class in joint models)
# Define the run as a test run or real analysis
test_run <- FALSE
fit_ctrl <- list(
  chains = 4,
  iter   = if (test_run) 400 else 4000,
  warmup = if (test_run) 200 else 2000,
  cores = 4, threads = 2,
  seed = 42, adapt_delta = 0.95, max_treedepth = 12
)
backend    <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
refit_mode <- "on_change"       # brms file_refit: "never", "on_change", "always"


######## Export settings #####################################################
save_plots   <- TRUE
plot_formats <- c("png", "jpeg")
plot_dpi     <- 300
plot_width   <- 7
plot_height  <- 5
save_draws   <- TRUE            # full posterior draws (rows = draws, cols = distances)
save_ppc     <- TRUE            # posterior predictive check plot per model


###############################################################################
######## Helpers: general ####################################################
###############################################################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}

as_flag <- function(x) as.character(x) %in% c("TRUE", "true", "1")

sub_dir <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

prep_model_data <- function(d) {
  d %>%
    mutate(distance_f = factor(distance, levels = dist_levels)) %>%
    droplevels()
}

model_family_obj <- function() {
  switch(model_family, gaussian = gaussian(), student = student(),
         stop("model_family must be 'gaussian' or 'student'"))
}


######## Helpers: self-documentation #########################################
settings_lines <- function() {
  c("## Settings",
    paste0("Input: ", input_file),
    paste0("Analysis tag: ", analysis_tag),
    paste0("Distance label: ", distance_col),
    paste0("Raw-value bounds applied: ", apply_bounds),
    paste0("Vessel subsample per subject x region: ",
           ifelse(is.null(vessel_subsample_n), "none (all vessels)", vessel_subsample_n),
           " (seed ", subsample_seed, ")"),
    paste0("Diameter mode: ", diameter_mode,
           if (diameter_mode == "matched")
             paste0(" (EPVS quantiles ", paste(matched_quantiles, collapse = "-"), ")") else ""),
    paste0("Family: ", model_family, " | heteroscedastic sigma: ", heteroscedastic),
    paste0("Backend: ", backend, " | refit mode: ", refit_mode),
    paste0("Sampling: ", paste(names(fit_ctrl), unlist(fit_ctrl), sep = " = ", collapse = ", ")),
    paste0("Low-support flag: < ", min_objects_flag, " objects per distance"))
}

doc_init <- function(dir, title, extra = character()) {
  f <- file.path(dir, "model_info.txt")
  writeLines(c(paste("#", title),
               paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               extra, "", settings_lines(), "",
               "## Session", capture.output(sessionInfo()), ""), f)
  f
}

doc_add <- function(f, ...) cat(c(...), file = f, sep = "\n", append = TRUE)

doc_capture <- function(f, x) doc_add(f, capture.output(print(x)))

doc_model <- function(doc, label, fit, data, dg, minutes) {
  doc_add(doc, "", paste0("### Model: ", label),
          paste0("Rows: ", nrow(data), " | objects: ", n_distinct(data$object_id),
                 " | subjects: ", paste(sort(unique(data$subjectID)), collapse = ", ")),
          sprintf("Fit/load time: %.1f min", minutes), "Formula:")
  doc_capture(doc, formula(fit))
  doc_add(doc, "Priors:")
  doc_capture(doc, prior_summary(fit))
  doc_add(doc, "Diagnostics:")
  doc_capture(doc, as.data.frame(dg))
}


######## Helpers: model specification ########################################
build_formula <- function(resp, data, joint) {
  rhs <- if (joint) "class * distance_f" else "distance_f"
  if (n_distinct(data$subjectID) >= 2) rhs <- paste(rhs, "+ subjectID")
  if (diameter_mode == "covariate")    rhs <- paste(rhs, "+ s_logdiam")
  re  <- if (joint) "(1 | gr(object_id, by = class))" else "(1 | object_id)"
  f_mu <- as.formula(paste(resp, "~", rhs, "+", re))
  if (!heteroscedastic) return(bf(f_mu))
  f_sigma <- as.formula(if (joint) "sigma ~ class + distance_f" else "sigma ~ distance_f")
  bf(f_mu, f_sigma)
}

# Weakly informative priors scaled to the response
build_priors <- function(data, resp) {
  s <- sd(data[[resp]], na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- 1
  p <- set_prior(sprintf("normal(0, %.4g)", 2.5 * s), class = "b") +
    set_prior(sprintf("student_t(3, 0, %.4g)", s), class = "sd")
  if (heteroscedastic) p <- p + set_prior("normal(0, 1)", class = "b", dpar = "sigma")
  p
}


######## Helpers: fitting and diagnostics ####################################
diag_log <- list()

model_diagnostics <- function(fit) {
  np  <- nuts_params(fit)
  div <- sum(np$Value[np$Parameter == "divergent__"])
  tdh <- sum(np$Value[np$Parameter == "treedepth__"] >= fit_ctrl$max_treedepth)
  sm  <- rstan::summary(fit$fit)$summary
  sm  <- sm[rownames(sm) != "lp__", , drop = FALSE]
  tibble(
    divergences    = div,
    treedepth_hits = tdh,
    max_rhat       = max(sm[, "Rhat"], na.rm = TRUE),
    worst_rhat_par = rownames(sm)[which.max(sm[, "Rhat"])],
    min_n_eff      = min(sm[, "n_eff"], na.rm = TRUE),
    worst_neff_par = rownames(sm)[which.min(sm[, "n_eff"])]
  )
}

fit_model <- function(data, resp, joint, file_stem, label, doc) {
  message(sprintf("  [%s] %d subject(s), %d objects, %d rows",
                  label, n_distinct(data$subjectID), n_distinct(data$object_id), nrow(data)))
  t0  <- Sys.time()
  fit <- brm(
    formula = build_formula(resp, data, joint),
    data    = data,
    family  = model_family_obj(),
    prior   = build_priors(data, resp),
    chains  = fit_ctrl$chains, iter = fit_ctrl$iter, warmup = fit_ctrl$warmup,
    cores   = fit_ctrl$cores, seed = fit_ctrl$seed,
    threads = if (backend == "cmdstanr" && fit_ctrl$threads > 1) threading(fit_ctrl$threads) else NULL,
    control = list(adapt_delta = fit_ctrl$adapt_delta, max_treedepth = fit_ctrl$max_treedepth),
    backend = backend,
    file = file_stem, file_refit = refit_mode,     # saves/reloads <file_stem>.rds
    silent = 2, refresh = 0
  )
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  
  dg <- model_diagnostics(fit)
  bad <- dg$divergences > 0 || dg$max_rhat > 1.01 || dg$min_n_eff < 400
  message(sprintf("    %.1f min | divergences %d | max Rhat %.3f (%s) | min n_eff %.0f (%s)%s",
                  minutes, dg$divergences, dg$max_rhat, dg$worst_rhat_par,
                  dg$min_n_eff, dg$worst_neff_par, if (bad) "  <-- CHECK" else ""))
  
  write_csv(dg, paste0(file_stem, "_diagnostics.csv"))
  diag_log[[length(diag_log) + 1]] <<- dg %>% mutate(model = label, file = file_stem, .before = 1)
  write_csv(bind_rows(diag_log), file.path(root_dir, "all_model_diagnostics.csv"))
  
  if (save_ppc) {
    p <- pp_check(fit, ndraws = 50) + ggtitle(label) + theme_minimal()
    save_plot(p, paste0(file_stem, "_ppc"))
  }
  doc_model(doc, label, fit, data, dg, minutes)
  fit
}


######## Helpers: posterior curves ###########################################
# Draws x distance matrix per class: expected value for a typical object
# (re_formula = NA), at the mean diameter (s_logdiam = 0), averaged equally
# over the subjects in the data.
epred_curves <- function(fit, data) {
  out <- list()
  for (cl in levels(data$class)) {
    d_present <- sort(unique(data$distance[data$class == cl]))
    nd <- expand.grid(distance = d_present, subjectID = unique(data$subjectID),
                      stringsAsFactors = FALSE)
    nd$distance_f <- factor(nd$distance, levels = levels(data$distance_f))
    nd$class      <- factor(cl, levels = levels(data$class))
    nd$s_logdiam  <- 0
    ep <- posterior_epred(fit, newdata = nd, re_formula = NA)
    m  <- sapply(d_present, function(dd) rowMeans(ep[, nd$distance == dd, drop = FALSE]))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(ep))
    colnames(m) <- as.character(d_present)
    out[[cl]] <- m
  }
  out
}

diff_matrix <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NULL)
  common <- intersect(colnames(a), colnames(b))
  common <- common[order(as.numeric(common))]
  if (length(common) == 0) return(NULL)
  n <- min(nrow(a), nrow(b))
  a[seq_len(n), common, drop = FALSE] - b[seq_len(n), common, drop = FALSE]
}

# Posterior summary per distance, with convergence checked on the reported
# quantity itself (draws are stored chain by chain)
summarise_draw_matrix <- function(m, n_chains, prob = 0.95) {
  a <- (1 - prob) / 2
  map_dfr(colnames(m), function(d) {
    v  <- m[, d]
    mm <- matrix(v, ncol = n_chains)
    tibble(distance = as.numeric(d), mean = mean(v), sd = sd(v),
           lower = unname(quantile(v, a)), upper = unname(quantile(v, 1 - a)),
           p_gt_0 = mean(v > 0),
           rhat = posterior::rhat(mm), ess_bulk = posterior::ess_bulk(mm))
  })
}

write_draws <- function(m, path) {
  if (!save_draws || is.null(m)) return(invisible(NULL))
  out <- as.data.frame(m)
  names(out) <- paste0("d", colnames(m))
  write_csv(cbind(draw = seq_len(nrow(out)), out), path)
}


######## Helpers: raw summaries and object counts ############################
raw_summary <- function(data, resp) {
  data %>%
    filter(!is.na(.data[[resp]])) %>%
    group_by(class, distance) %>%
    summarise(mean = mean(.data[[resp]]), sd = sd(.data[[resp]]),
              n_objects = n_distinct(object_id), n_subjects = n_distinct(subjectID),
              .groups = "drop") %>%
    mutate(class = as.character(class),
           se = sd / sqrt(n_objects),
           low_support = n_objects < min_objects_flag)
}

count_objects <- function(data, resp) {
  data %>%
    filter(!is.na(.data[[resp]])) %>%
    count(Stage, subjectID, class, distance, name = "n_objects") %>%
    arrange(class, subjectID, distance)
}

raw_diff <- function(r_k, r_0) {
  inner_join(r_0 %>% select(distance, mean_0 = mean, se_0 = se, nobj_0 = n_objects),
             r_k %>% select(distance, mean_k = mean, se_k = se, nobj_k = n_objects),
             by = "distance") %>%
    mutate(diff = mean_k - mean_0, se_diff = sqrt(se_k^2 + se_0^2))
}

raw_epvs_minus_vessel <- function(r) {
  e <- r %>% filter(class == "EPVS")
  v <- r %>% filter(class == "vessel")
  if (nrow(e) == 0 || nrow(v) == 0) return(NULL)
  raw_diff(e, v)
}


######## Helpers: plots ######################################################
save_plot <- function(p, stem_path, width = plot_width, height = plot_height) {
  if (!save_plots) return(invisible(NULL))
  for (fmt in plot_formats) {
    ggsave(paste0(stem_path, ".", fmt), p, width = width, height = height,
           dpi = plot_dpi, bg = "white")
  }
}

x_lab <- "Distance from boundary (\u00b5m, shell outer edge)"

plot_curves <- function(fit_df, raw_df, color_col, colors, title, ylab,
                        facet_col = NULL, zero_line = FALSE) {
  p <- ggplot()
  if (zero_line) p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
  p <- p +
    geom_ribbon(data = fit_df, aes(x = distance, ymin = lower, ymax = upper,
                                   fill = .data[[color_col]]), alpha = 0.18) +
    geom_line(data = fit_df, aes(x = distance, y = mean, color = .data[[color_col]]),
              linewidth = 0.8)
  if (!is.null(raw_df) && nrow(raw_df) > 0) {
    raw_df$low_support <- factor(raw_df$low_support, levels = c(FALSE, TRUE))
    p <- p +
      geom_errorbar(data = raw_df, aes(x = distance, ymin = mean - se, ymax = mean + se,
                                       color = .data[[color_col]]), width = 8, alpha = 0.5) +
      geom_point(data = raw_df, aes(x = distance, y = mean, color = .data[[color_col]],
                                    shape = low_support), size = 1.6) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), drop = FALSE, name = NULL,
                         labels = c(`FALSE` = paste0("\u2265 ", min_objects_flag, " objects"),
                                    `TRUE`  = paste0("< ",      min_objects_flag, " objects")))
  }
  p <- p + (if (is.null(colors))
    scale_color_discrete(aesthetics = c("color", "fill"), name = NULL)
    else
      scale_color_manual(values = colors, aesthetics = c("color", "fill"), name = NULL))
  if (!is.null(facet_col)) p <- p + facet_wrap(facet_col, nrow = 1)
  p + scale_x_continuous(breaks = dist_levels) +
    labs(x = x_lab, y = ylab, title = title,
         subtitle = "lines/ribbons = posterior mean (95% CrI); points = raw mean across objects (\u00b1 SE)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_diff <- function(model_df, raw_df, title, ylab, facet_col = NULL) {
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_ribbon(data = model_df, aes(x = distance, ymin = lower, ymax = upper),
                alpha = 0.2, fill = "steelblue") +
    geom_line(data = model_df, aes(x = distance, y = mean), color = "steelblue", linewidth = 0.8)
  if (!is.null(raw_df) && nrow(raw_df) > 0) {
    p <- p +
      geom_errorbar(data = raw_df, aes(x = distance, ymin = diff - se_diff, ymax = diff + se_diff),
                    color = "firebrick", width = 8, alpha = 0.7) +
      geom_point(data = raw_df, aes(x = distance, y = diff), color = "firebrick", size = 1.3)
  }
  if (!is.null(facet_col)) p <- p + facet_wrap(facet_col, nrow = 1)
  p + scale_x_continuous(breaks = dist_levels) +
    labs(x = x_lab, y = ylab, title = title,
         subtitle = "blue = model difference (95% CrI); red = raw-mean difference (\u00b1 SE)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1))
}


###############################################################################
######## Analysis blocks #####################################################
###############################################################################

# A) and B): one model per stage, stage contrasts, combined outputs.
# joint = FALSE : data contain one class (EPVS or vessel)
# joint = TRUE  : data contain both classes
run_stage_set <- function(dat, resp, joint, block, stem, doc, ylab, label_prefix) {
  stages_here <- intersect(stage_levels, unique(as.character(dat$Stage)))
  if (!"control" %in% stages_here) {
    message("  No control data - skipping ", label_prefix)
    doc_add(doc, "", "No control data: block skipped.")
    return(invisible(NULL))
  }
  
  res <- list()
  for (s in stages_here) {
    ds <- dat %>% filter(Stage == s) %>% prep_model_data()
    if (joint && nlevels(ds$class) < 2) {
      message("  [", label_prefix, "/", s, "] only one class present - joint model skipped")
      doc_add(doc, "", paste0("Stage ", s, ": only one class present, joint model skipped."))
      next
    }
    sdir <- sub_dir(block, s)
    sst  <- file.path(sdir, paste0(stem, "_", s))
    fit  <- fit_model(ds, resp, joint, paste0(sst, "_model"), paste0(label_prefix, "/", s), doc)
    nc   <- nchains(fit)
    
    curves  <- epred_curves(fit, ds)
    raw_c   <- raw_summary(ds, resp)
    fit_tbl <- imap_dfr(curves, ~ summarise_draw_matrix(.x, nc) %>% mutate(class = .y, .before = 1)) %>%
      left_join(raw_c %>% select(class, distance, n_objects, n_subjects, low_support),
                by = c("class", "distance"))
    
    write_csv(fit_tbl, paste0(sst, "_fitted.csv"))
    write_csv(raw_c,   paste0(sst, "_raw.csv"))
    write_csv(count_objects(ds, resp), paste0(sst, "_objects_per_distance.csv"))
    iwalk(curves, ~ write_draws(.x, paste0(sst, "_", .y, "_draws.csv")))
    
    ev <- NULL; ev_tbl <- NULL; ev_raw <- NULL
    if (joint) {
      ev     <- diff_matrix(curves[["EPVS"]], curves[["vessel"]])
      ev_tbl <- summarise_draw_matrix(ev, nc)
      ev_raw <- raw_epvs_minus_vessel(raw_c)
      write_csv(ev_tbl, paste0(sst, "_EPVS_minus_vessel.csv"))
      if (!is.null(ev_raw)) write_csv(ev_raw, paste0(sst, "_EPVS_minus_vessel_raw.csv"))
      write_draws(ev, paste0(sst, "_EPVS_minus_vessel_draws.csv"))
      
      p <- plot_curves(fit_tbl, raw_c, "class", class_colors,
                       paste0(label_prefix, " / ", s, ": EPVS and vessel"), ylab)
      save_plot(p, paste0(sst, "_curves"))
    }
    res[[s]] <- list(curves = curves, raw = raw_c, fit_tbl = fit_tbl,
                     ev = ev, ev_tbl = ev_tbl, ev_raw = ev_raw, nc = nc)
  }
  if (is.null(res$control)) return(invisible(NULL))
  
  # ---- Stage contrasts vs control --------------------------------------------
  diff_all <- list(); rawdiff_all <- list(); did_all <- list()
  for (k in intersect(nonzero_stages, names(res))) {
    kdir <- sub_dir(block, paste0(k, "_vs_control"))
    kst  <- file.path(kdir, paste0(stem, "_", k, "_vs_control"))
    nc   <- res[[k]]$nc
    
    for (cl in intersect(names(res[[k]]$curves), names(res$control$curves))) {
      dm <- diff_matrix(res[[k]]$curves[[cl]], res$control$curves[[cl]])
      if (is.null(dm)) next
      dt <- summarise_draw_matrix(dm, nc) %>% mutate(class = cl, contrast = paste(k, "- control"), .before = 1)
      rd <- raw_diff(res[[k]]$raw %>% filter(class == cl), res$control$raw %>% filter(class == cl)) %>%
        mutate(class = cl, contrast = paste(k, "- control"), .before = 1)
      write_csv(dt, paste0(kst, "_", cl, "_model_diff.csv"))
      write_csv(rd, paste0(kst, "_", cl, "_raw_diff.csv"))
      write_draws(dm, paste0(kst, "_", cl, "_diff_draws.csv"))
      diff_all[[paste(k, cl)]]    <- dt
      rawdiff_all[[paste(k, cl)]] <- rd
    }
    
    if (joint && !is.null(res[[k]]$ev) && !is.null(res$control$ev)) {
      did <- diff_matrix(res[[k]]$ev, res$control$ev)
      if (!is.null(did)) {
        dd <- summarise_draw_matrix(did, nc) %>%
          mutate(contrast = paste0("(EPVS - vessel) ", k, " - control"), .before = 1)
        write_csv(dd, paste0(kst, "_EPVS_minus_vessel_diff.csv"))
        write_draws(did, paste0(kst, "_EPVS_minus_vessel_diff_draws.csv"))
        did_all[[k]] <- dd
      }
    }
  }
  
  # ---- Combined outputs ------------------------------------------------------
  adir <- sub_dir(block, "all_stages")
  ast  <- file.path(adir, paste0(stem, "_all_stages"))
  fit_all <- map_dfr(res, "fit_tbl", .id = "stage") %>% mutate(stage = factor(stage, stage_levels))
  raw_all <- map_dfr(res, "raw",     .id = "stage") %>% mutate(stage = factor(stage, stage_levels))
  write_csv(fit_all, paste0(ast, "_fitted.csv"))
  write_csv(raw_all, paste0(ast, "_raw.csv"))
  
  if (joint) {
    save_plot(plot_curves(fit_all, raw_all, "class", class_colors,
                          paste0(label_prefix, ": EPVS and vessel by stage"), ylab, facet_col = "stage"),
              paste0(ast, "_curves"), width = max(plot_width, plot_width * 0.6 * length(res)))
    
    ev_all  <- map_dfr(res, "ev_tbl", .id = "stage") %>% mutate(stage = factor(stage, stage_levels))
    evr_all <- map_dfr(res, "ev_raw", .id = "stage")
    if (nrow(evr_all) > 0) {
      evr_all <- evr_all %>%
        transmute(stage = factor(stage, stage_levels), distance, mean = diff, se = se_diff,
                  low_support = pmin(nobj_0, nobj_k) < min_objects_flag)
    }
    write_csv(ev_all, paste0(ast, "_EPVS_minus_vessel.csv"))
    save_plot(plot_curves(ev_all, evr_all, "stage", stage_colors,
                          paste0(label_prefix, ": EPVS - vessel by stage"),
                          paste0("EPVS - vessel, ", ylab), zero_line = TRUE),
              paste0(ast, "_EPVS_minus_vessel"))
  } else {
    save_plot(plot_curves(fit_all, raw_all, "stage", stage_colors,
                          paste0(label_prefix, ": fitted curve by stage"), ylab),
              paste0(ast, "_curves"))
  }
  
  if (length(diff_all) > 0) {
    d_all  <- bind_rows(diff_all)
    rd_all <- bind_rows(rawdiff_all)
    write_csv(d_all,  paste0(ast, "_model_diff.csv"))
    write_csv(rd_all, paste0(ast, "_raw_diff.csv"))
    d_all$panel  <- paste(d_all$class,  d_all$contrast,  sep = ": ")
    rd_all$panel <- paste(rd_all$class, rd_all$contrast, sep = ": ")
    save_plot(plot_diff(d_all, rd_all, paste0(label_prefix, ": stage - control"),
                        paste0("Difference, ", ylab), facet_col = "panel"),
              paste0(ast, "_stage_diffs"),
              width = max(plot_width, plot_width * 0.6 * max(1, n_distinct(d_all$panel))))
  }
  if (length(did_all) > 0) {
    did_df <- bind_rows(did_all)
    write_csv(did_df, paste0(ast, "_EPVS_minus_vessel_diff.csv"))
    save_plot(plot_diff(did_df, NULL, paste0(label_prefix, ": change in EPVS - vessel vs control"),
                        paste0("Difference in (EPVS - vessel), ", ylab), facet_col = "contrast"),
              paste0(ast, "_EPVS_minus_vessel_diff"),
              width = max(plot_width, plot_width * 0.6 * max(1, length(did_all))))
  }
  invisible(res)
}

# C) one model per subject (joint if both classes are present)
run_subject_block <- function(ds_all, resp, block, stem, doc, ylab, label) {
  ds <- prep_model_data(ds_all)
  joint <- nlevels(ds$class) == 2
  sst <- file.path(block, stem)
  fit <- fit_model(ds, resp, joint, paste0(sst, "_model"), label, doc)
  nc  <- nchains(fit)
  
  curves  <- epred_curves(fit, ds)
  raw_c   <- raw_summary(ds, resp)
  fit_tbl <- imap_dfr(curves, ~ summarise_draw_matrix(.x, nc) %>% mutate(class = .y, .before = 1)) %>%
    left_join(raw_c %>% select(class, distance, n_objects, low_support), by = c("class", "distance"))
  write_csv(fit_tbl, paste0(sst, "_fitted.csv"))
  write_csv(raw_c,   paste0(sst, "_raw.csv"))
  write_csv(count_objects(ds, resp), paste0(sst, "_objects_per_distance.csv"))
  iwalk(curves, ~ write_draws(.x, paste0(sst, "_", .y, "_draws.csv")))
  save_plot(plot_curves(fit_tbl, raw_c, "class", class_colors, paste0(label, ": EPVS and vessel"), ylab),
            paste0(sst, "_curves"))
  
  ev_tbl <- NULL
  if (joint) {
    ev     <- diff_matrix(curves[["EPVS"]], curves[["vessel"]])
    ev_tbl <- summarise_draw_matrix(ev, nc)
    write_csv(ev_tbl, paste0(sst, "_EPVS_minus_vessel.csv"))
    write_draws(ev, paste0(sst, "_EPVS_minus_vessel_draws.csv"))
  }
  list(fit_tbl = fit_tbl, ev_tbl = ev_tbl)
}


###############################################################################
######## BEGIN LOGIC #########################################################
###############################################################################
input_file   <- file.path(common_dir, input_filename)
base_name    <- tools::file_path_sans_ext(input_filename)
analysis_tag <- paste0(if (diameter_mode == "none") "main" else paste0("diam_", diameter_mode),
                       if (is.null(vessel_subsample_n)) "" else paste0("_vsub", vessel_subsample_n),
                       if (test_run) "_TEST" else "")
root_dir     <- sub_dir(output_filepath, base_name, "brms", analysis_tag)
message("ANALYZING: ", input_filename, "\nOUTPUT:    ", root_dir, "\nBACKEND:   ", backend)

######## 1) Import ###########################################################
raw_in <- read_csv(input_file, show_col_types = FALSE)
needed <- c("subject", "severity", "region", "class", "object_id", distance_col,
            "enough_voxels", "obj_diam_um", "mus_ok", "ret_ok", unlist(response_cols))
missing_cols <- setdiff(needed, names(raw_in))
if (length(missing_cols) > 0) stop("CSV is missing columns: ", paste(missing_cols, collapse = ", "))

keep_row <- as_flag(raw_in$enough_voxels)
message("Rows with enough voxels: ", sum(keep_row), " of ", nrow(raw_in))

shells <- raw_in[keep_row, ] %>%
  transmute(
    subjectID = as.character(subject),
    Stage     = factor(severity, levels = stage_levels),
    Region    = normalize_region(region),
    class     = factor(as.character(class), levels = class_levels),
    object_id = as.character(object_id),
    distance  = .data[[distance_col]],
    obj_diam_um, n_vox,
    mus_ok = as_flag(mus_ok), ret_ok = as_flag(ret_ok),
    scattering_raw = mus, scattering_pct = mus_pct,
    retardance_raw = ret, retardance_pct = ret_pct
  )
n_bad <- sum(is.na(shells$Stage) | is.na(shells$class))
if (n_bad > 0) warning(n_bad, " rows with unexpected severity or class labels dropped")
shells <- shells %>% filter(!is.na(Stage), !is.na(class))

if (apply_bounds) {
  message("Rows failing bounds: scattering ", sum(!shells$mus_ok),
          ", retardance ", sum(!shells$ret_ok))
  shells <- shells %>%
    mutate(across(c(scattering_raw, scattering_pct), ~ if_else(mus_ok, .x, NA_real_)),
           across(c(retardance_raw, retardance_pct), ~ if_else(ret_ok, .x, NA_real_)))
}

######## 2) Diameter sensitivity #############################################
obj_tbl <- shells %>% distinct(object_id, Region, subjectID, class, obj_diam_um)

if (diameter_mode == "matched") {
  rng <- obj_tbl %>% filter(class == "EPVS") %>% group_by(Region) %>%
    summarise(lo = quantile(obj_diam_um, matched_quantiles[1]),
              hi = quantile(obj_diam_um, matched_quantiles[2]), .groups = "drop")
  keep_ids <- obj_tbl %>% left_join(rng, by = "Region") %>%
    filter(class == "EPVS" | is.na(lo) | (obj_diam_um >= lo & obj_diam_um <= hi)) %>%
    pull(object_id)
  message("Diameter matching: kept ", sum(obj_tbl$class == "vessel" & obj_tbl$object_id %in% keep_ids),
          " of ", sum(obj_tbl$class == "vessel"), " vessel objects")
  write_csv(rng, file.path(root_dir, "diameter_matching_ranges.csv"))
  shells  <- shells  %>% filter(object_id %in% keep_ids)
  obj_tbl <- obj_tbl %>% filter(object_id %in% keep_ids)
}

if (diameter_mode == "covariate") {
  ld <- log(obj_tbl$obj_diam_um)
  shells$s_logdiam <- (log(shells$obj_diam_um) - mean(ld)) / sd(ld)
} else {
  shells$s_logdiam <- 0
}

######## 3) Vessel subsampling ###############################################
if (!is.null(vessel_subsample_n)) {
  set.seed(subsample_seed)
  ves_ids <- shells %>% filter(class == "vessel") %>%
    distinct(subjectID, Region, object_id) %>%
    arrange(subjectID, Region, object_id) %>%
    group_by(subjectID, Region) %>%
    slice_sample(n = vessel_subsample_n) %>%
    ungroup()
  write_csv(ves_ids, file.path(root_dir, "vessel_subsample_object_ids.csv"))
  shells <- shells %>% filter(class == "EPVS" | object_id %in% ves_ids$object_id)
}

dist_levels <- sort(unique(shells$distance))
message("Distances: ", paste(dist_levels, collapse = ", "))

obj_counts <- shells %>% group_by(Stage, Region, class, subjectID) %>%
  summarise(objects = n_distinct(object_id), rows = n(), .groups = "drop")
print(obj_counts, n = Inf)
write_csv(obj_counts, file.path(root_dir, "objects_per_subject_region_class.csv"))
writeLines(c("# Run settings", paste("Generated:", format(Sys.time())), "", settings_lines()),
           file.path(root_dir, "run_settings.txt"))

######## 4) Main loop ########################################################
for (scl in scales) {
  for (prop in names(response_cols)) {
    resp  <- paste0(prop, "_", scl)
    y_lab <- scale_ylab[[scl]][[prop]]
    
    for (rgn in regions) {
      reg_dat <- shells %>% filter(Region == rgn, !is.na(.data[[resp]]))
      if (nrow(reg_dat) == 0) next
      reg_dir <- file.path(root_dir, scl, prop, rgn)
      tag     <- paste(scl, prop, rgn, sep = "/")
      message("\n===== ", tag, " =====")
      
      ## A) Group level, separate models ---------------------------------------
      if (run_separate) {
        for (src in names(sources)) {
          d_src <- reg_dat %>% filter(class == sources[[src]])
          if (nrow(d_src) == 0) next
          block <- sub_dir(reg_dir, "group", src)
          doc   <- doc_init(block, paste0("Group-level separate model: ", tag, "/", src),
                            "Model per stage: response ~ distance_f [+ subjectID] + (1 | object_id)")
          run_stage_set(d_src, resp, joint = FALSE, block,
                        stem = paste(scl, prop, rgn, src, sep = "_"), doc, y_lab,
                        label_prefix = paste0(tag, "/group/", src))
        }
      }
      
      ## B) Group level, joint EPVS + vessel ------------------------------------
      if (run_joint) {
        block <- sub_dir(reg_dir, "group", "joint")
        doc   <- doc_init(block, paste0("Group-level joint EPVS + vessel model: ", tag),
                          c("Model per stage: response ~ class * distance_f [+ subjectID] + (1 | gr(object_id, by = class))",
                            "EPVS - vessel = difference of class-specific posterior curves, draw by draw"))
        run_stage_set(reg_dat, resp, joint = TRUE, block,
                      stem = paste(scl, prop, rgn, "joint", sep = "_"), doc, y_lab,
                      label_prefix = paste0(tag, "/group/joint"))
      }
      
      ## C) Subject level --------------------------------------------------------
      if (run_subject) {
        subj_fit <- list(); subj_ev <- list()
        for (subj in sort(unique(reg_dat$subjectID))) {
          d_sub <- reg_dat %>% filter(subjectID == subj)
          stg   <- as.character(unique(d_sub$Stage))
          block <- sub_dir(reg_dir, "subject", subj)
          doc   <- doc_init(block, paste0("Subject-level model: ", tag, "/", subj, " (", stg, ")"),
                            "Model: response ~ class * distance_f + (1 | gr(object_id, by = class)), one subject")
          out <- run_subject_block(d_sub, resp, block, paste(scl, prop, rgn, subj, sep = "_"),
                                   doc, y_lab, paste0(tag, "/subject/", subj))
          subj_fit[[subj]] <- out$fit_tbl %>% mutate(subjectID = subj, stage = stg, .before = 1)
          if (!is.null(out$ev_tbl)) subj_ev[[subj]] <- out$ev_tbl %>% mutate(subjectID = subj, stage = stg, .before = 1)
        }
        
        sdir <- sub_dir(reg_dir, "subject", "all_subjects")
        sst  <- file.path(sdir, paste(scl, prop, rgn, "all_subjects", sep = "_"))
        sf   <- bind_rows(subj_fit)
        write_csv(sf, paste0(sst, "_fitted.csv"))
        sf$panel <- paste0(sf$subjectID, " (", sf$stage, ")")
        save_plot(plot_curves(sf, NULL, "class", class_colors,
                              paste0(tag, ": EPVS and vessel by subject"), y_lab, facet_col = "panel"),
                  paste0(sst, "_curves"), width = max(plot_width, plot_width * 0.55 * n_distinct(sf$panel)))
        if (length(subj_ev) > 0) {
          se <- bind_rows(subj_ev)
          write_csv(se, paste0(sst, "_EPVS_minus_vessel.csv"))
          save_plot(plot_curves(se, NULL, "subjectID", NULL,
                                paste0(tag, ": EPVS - vessel by subject"),
                                paste0("EPVS - vessel, ", y_lab), facet_col = "stage", zero_line = TRUE),
                    paste0(sst, "_EPVS_minus_vessel"),
                    width = max(plot_width, plot_width * 0.6 * n_distinct(se$stage)))
        }
      }
      message("  Finished: ", reg_dir)
    }
  }
}

message("\nAll done. Diagnostics for every model: ", file.path(root_dir, "all_model_diagnostics.csv"))
