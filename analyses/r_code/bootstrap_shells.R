###############################################################################
# Hierarchical bootstrap summaries of shell measurements (OCT or histology)
#
# A descriptive alternative to the brms models:
#   1. Object level   : each object's value in each 40 µm shell (from the CSV)
#   2. Subject level  : mean over that subject's objects, per class and
#                       distance; 95% interval from resampling OBJECTS within
#                       the subject (B bootstrap resamples)
#   3. Group level    : average of the subject curves (each subject counts
#                       equally); its bootstrap draws are the average of the
#                       subjects' bootstrap draws
#   4. Differences    : EPVS - vessel within each subject and stage,
#                       stage - control, and the change in EPVS - vessel,
#                       all computed draw by draw
#
# Intervals reflect sampling of objects within the measured subjects only;
# they do not include variation between subjects (n = 1-2 per stage).
#
# Output uses the same folder and file names as the brms scripts, under
# .../<base_name>/bootstrap/B<B>/, so plot_fig3_shells.m, plot_histo_shells.m
# and plot_joint_epvs_vessel.m can plot it with cfg.method = 'bootstrap'.
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

######## Settings ############################################################
dataset <- "histo"                 # "oct" or "histo"

B                <- 2000         # bootstrap resamples
boot_seed        <- 2026
ci_prob          <- 0.95
min_objects_flag <- 10           # flag distances with fewer objects
distance_col     <- "d_outer_um" # label shells by outer edge (40 ... 520)
stage_levels     <- c("control", "moderate", "severe")
nonzero_stages   <- c("moderate", "severe")
regions          <- c("front", "occip")
save_draws       <- TRUE         # bootstrap draws (rows = resamples, cols = distances)

if (dataset == "oct") {
  common_dir      <- "/projectnb/npbssmic/ns/CAA"
  input_filename  <- "parenchyma_shells_2026-10-07_pct.csv"
  output_filepath <- "/projectnb/npbssmic/ns/CAA/beta_stats"
  scales          <- c("pct")                     # "raw" and/or "pct"
  props           <- c("scattering", "retardance")
  apply_bounds    <- TRUE                         # mus_ok / ret_ok flags
} else if (dataset == "histo") {
  common_dir      <- "/projectnb/npbssmic/ns/CAA/histology"
  input_filename  <- "histo_shells_2026-10-08.csv"
  output_filepath <- "/projectnb/npbssmic/ns/CAA/histology/beta_stats"
  scales          <- c("z")                       # "z", "hmatch" and/or "raw"
  props           <- c("LHE", "CD68", "GFAP")     # stains
} else {
  stop("dataset must be 'oct' or 'histo'")
}


######## Helpers #############################################################
normalize_region <- function(r) {
  r2 <- tolower(as.character(r))
  ifelse(str_detect(r2, "occip"), "occip",
         ifelse(str_detect(r2, "front|frontal"), "front", as.character(r)))
}
as_flag <- function(x) as.character(x) %in% c("TRUE", "true", "1")
sub_dir <- function(...) { d <- file.path(...); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

# Bootstrap the per-distance mean over one subject's objects of one class.
# Objects are resampled with replacement as whole rows (all their shells).
boot_curve <- function(d) {
  w <- d %>% select(object_id, distance, value) %>%
    pivot_wider(names_from = distance, values_from = value, values_fn = mean)
  Y <- matrix(NA_real_, nrow(w), length(dist_levels),
              dimnames = list(NULL, as.character(dist_levels)))
  cols <- intersect(colnames(Y), names(w))
  Y[, cols] <- as.matrix(w[, cols])
  M  <- !is.na(Y) * 1
  Y0 <- Y; Y0[is.na(Y0)] <- 0
  n   <- nrow(Y)
  cnt <- colSums(M)
  est <- ifelse(cnt > 0, colSums(Y0) / cnt, NA_real_)
  W   <- t(rmultinom(B, n, rep(1, n)))         # B x n resampling counts
  bm  <- (W %*% Y0) / (W %*% M)                # B x D bootstrap means
  bm[!is.finite(bm)] <- NA
  colnames(bm) <- colnames(Y)
  list(est = est, boot = bm, n_objects = cnt, n_subjects = as.numeric(cnt > 0))
}

# Equal-weight average of several subjects' curves (draw by draw)
avg_curves <- function(lst) {
  est_m <- do.call(cbind, map(lst, "est"))
  n_sub <- rowSums(!is.na(est_m))
  est   <- ifelse(n_sub > 0, rowSums(est_m, na.rm = TRUE) / n_sub, NA_real_)
  sb <- Reduce(`+`, map(lst, ~ replace(.x$boot, is.na(.x$boot), 0)))
  cb <- Reduce(`+`, map(lst, ~ (!is.na(.x$boot)) * 1))
  bm <- sb / cb; bm[!is.finite(bm)] <- NA
  list(est = est, boot = bm, n_objects = Reduce(`+`, map(lst, "n_objects")), n_subjects = n_sub)
}

diff_curves <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NULL)
  list(est = a$est - b$est, boot = a$boot - b$boot,
       n_objects = pmin(a$n_objects, b$n_objects),
       n_subjects = pmin(a$n_subjects, b$n_subjects))
}

summ <- function(x) {
  a <- (1 - ci_prob) / 2
  tibble(distance   = dist_levels,
         mean       = unname(x$est),
         lower      = apply(x$boot, 2, quantile, probs = a,     na.rm = TRUE, names = FALSE),
         upper      = apply(x$boot, 2, quantile, probs = 1 - a, na.rm = TRUE, names = FALSE),
         sd_boot    = apply(x$boot, 2, sd, na.rm = TRUE),
         frac_gt_0  = colMeans(x$boot > 0, na.rm = TRUE),
         n_objects  = unname(x$n_objects),
         n_subjects = unname(x$n_subjects)) %>%
    filter(!is.na(mean)) %>%
    mutate(low_support = n_objects < min_objects_flag)
}

write_draws <- function(x, path) {
  if (!save_draws || is.null(x)) return(invisible(NULL))
  out <- as.data.frame(x$boot)
  names(out) <- paste0("d", colnames(x$boot))
  write_csv(cbind(draw = seq_len(nrow(out)), out), path)
}


######## Import ##############################################################
input_file   <- file.path(common_dir, input_filename)
base_name    <- tools::file_path_sans_ext(input_filename)
analysis_tag <- paste0("B", B)
root_dir     <- sub_dir(output_filepath, base_name, "bootstrap", analysis_tag)
message("DATASET: ", dataset, "\nINPUT:   ", input_file, "\nOUTPUT:  ", root_dir)

csv <- read_csv(input_file, show_col_types = FALSE)

if (dataset == "oct") {
  long <- csv[as_flag(csv$enough_voxels), ] %>%
    transmute(subjectID = as.character(subject),
              Stage     = factor(severity, levels = stage_levels),
              Region    = normalize_region(region),
              class     = as.character(class),
              object_id = as.character(object_id),
              distance  = .data[[distance_col]],
              mus_ok = as_flag(mus_ok), ret_ok = as_flag(ret_ok),
              scattering_raw = mus, scattering_pct = mus_pct,
              retardance_raw = ret, retardance_pct = ret_pct) %>%
    pivot_longer(c(scattering_raw, scattering_pct, retardance_raw, retardance_pct),
                 names_to = c("prop", "scale"), names_sep = "_", values_to = "value")
  if (apply_bounds) {
    long <- long %>% filter(!(prop == "scattering" & !mus_ok), !(prop == "retardance" & !ret_ok))
  }
  long <- long %>% select(-mus_ok, -ret_ok)
} else {
  long <- csv[as_flag(csv$enough_pixels), ] %>%
    transmute(subjectID = as.character(subject),
              Stage     = factor(severity, levels = stage_levels),
              Region    = normalize_region(region),
              class     = as.character(class),
              object_id = as.character(object_id),
              distance  = .data[[distance_col]],
              prop      = as.character(stain),
              across(any_of(c("z", "hmatch", "raw")))) %>%
    pivot_longer(any_of(c("z", "hmatch", "raw")), names_to = "scale", values_to = "value")
}

long <- long %>%
  filter(!is.na(Stage), class %in% c("EPVS", "vessel"), Region %in% regions,
         scale %in% scales, prop %in% props, is.finite(value))
dist_levels <- sort(unique(long$distance))

counts <- long %>% distinct(scale, prop, Region, Stage, subjectID, class, object_id) %>%
  count(scale, prop, Region, Stage, subjectID, class, name = "objects")
print(counts, n = Inf)
write_csv(counts, file.path(root_dir, "objects_per_subject.csv"))

writeLines(c("# Hierarchical bootstrap summaries",
             paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             paste("Dataset:", dataset), paste("Input:", input_file),
             paste("Scales:", paste(scales, collapse = ", ")),
             paste("Properties/stains:", paste(props, collapse = ", ")),
             paste("Bootstrap resamples:", B, "| seed:", boot_seed, "| interval:", ci_prob),
             paste("Distance label:", distance_col, "| low-support flag: <", min_objects_flag, "objects"),
             if (dataset == "oct") paste("Raw-value bounds applied:", apply_bounds) else "",
             "",
             "Subject curve  = mean over the subject's objects at each distance;",
             "                 interval from resampling objects within the subject.",
             "Group curve    = equal-weight average of subject curves in a stage.",
             "Differences    = computed draw by draw (EPVS - vessel, stage - control,",
             "                 change in EPVS - vessel).",
             "Intervals do NOT include between-subject variation.",
             "File names mirror the brms outputs ('model_diff' = bootstrap difference).",
             "", "## Session", capture.output(sessionInfo())),
           file.path(root_dir, "method_info.txt"))


######## Main loop ###########################################################
for (scl in scales) {
  for (prp in props) {
    for (rgn in regions) {
      d_all <- long %>% filter(scale == scl, prop == prp, Region == rgn)
      if (nrow(d_all) == 0) next
      tag <- paste(scl, prp, rgn, sep = "/")
      message("\n===== ", tag, " =====")
      set.seed(boot_seed)
      reg_dir <- file.path(root_dir, scl, prp, rgn)

      subj_tbl <- d_all %>% distinct(subjectID, Stage) %>% arrange(Stage, subjectID)

      ## Subject-level bootstraps ---------------------------------------------
      sb <- list()
      for (s in subj_tbl$subjectID) {
        sb[[s]] <- list()
        for (cl in c("EPVS", "vessel")) {
          ds <- d_all %>% filter(subjectID == s, class == cl)
          if (nrow(ds) > 0) sb[[s]][[cl]] <- boot_curve(ds)
        }
      }

      subj_fit <- list(); subj_ev <- list()
      for (i in seq_len(nrow(subj_tbl))) {
        s <- subj_tbl$subjectID[i]; stg <- as.character(subj_tbl$Stage[i])
        for (cl in names(sb[[s]])) {
          subj_fit[[length(subj_fit) + 1]] <- summ(sb[[s]][[cl]]) %>%
            mutate(subjectID = s, stage = stg, class = cl, .before = 1)
        }
        ev <- diff_curves(sb[[s]][["EPVS"]], sb[[s]][["vessel"]])
        if (!is.null(ev)) subj_ev[[length(subj_ev) + 1]] <- summ(ev) %>%
          mutate(subjectID = s, stage = stg, .before = 1)
      }
      sdir <- sub_dir(reg_dir, "subject", "all_subjects")
      sst  <- file.path(sdir, paste(scl, prp, rgn, "all_subjects", sep = "_"))
      write_csv(bind_rows(subj_fit), paste0(sst, "_fitted.csv"))
      if (length(subj_ev) > 0) write_csv(bind_rows(subj_ev), paste0(sst, "_EPVS_minus_vessel.csv"))

      ## Group level: equal-weight average of subjects -------------------------
      grp <- list(); grp_ev <- list()
      for (stg in stage_levels) {
        ids <- subj_tbl$subjectID[subj_tbl$Stage == stg]
        if (length(ids) == 0) next
        grp[[stg]] <- list()
        for (cl in c("EPVS", "vessel")) {
          lst <- compact(map(ids, ~ sb[[.x]][[cl]]))
          if (length(lst) > 0) grp[[stg]][[cl]] <- avg_curves(lst)
        }
        grp_ev[[stg]] <- diff_curves(grp[[stg]][["EPVS"]], grp[[stg]][["vessel"]])
      }

      fit_tbl <- imap_dfr(grp, function(g, stg)
        imap_dfr(g, ~ summ(.x) %>% mutate(stage = stg, class = .y, .before = 1)))

      diff_tbl <- list(); did_tbl <- list()
      ddir <- sub_dir(reg_dir, "group", "joint", "draws")
      dst  <- file.path(ddir, paste(scl, prp, rgn, "joint", sep = "_"))
      for (stg in names(grp)) {
        iwalk(grp[[stg]], ~ write_draws(.x, paste0(dst, "_", stg, "_", .y, "_draws.csv")))
        write_draws(grp_ev[[stg]], paste0(dst, "_", stg, "_EPVS_minus_vessel_draws.csv"))
      }
      for (k in intersect(nonzero_stages, names(grp))) {
        if (is.null(grp[["control"]])) break
        for (cl in c("EPVS", "vessel")) {
          dk <- diff_curves(grp[[k]][[cl]], grp[["control"]][[cl]])
          if (is.null(dk)) next
          diff_tbl[[paste(k, cl)]] <- summ(dk) %>%
            mutate(class = cl, contrast = paste(k, "- control"), .before = 1)
          write_draws(dk, paste0(dst, "_", k, "_vs_control_", cl, "_diff_draws.csv"))
        }
        did <- diff_curves(grp_ev[[k]], grp_ev[["control"]])
        if (!is.null(did)) {
          did_tbl[[k]] <- summ(did) %>%
            mutate(contrast = paste0("(EPVS - vessel) ", k, " - control"), .before = 1)
          write_draws(did, paste0(dst, "_", k, "_vs_control_EPVS_minus_vessel_diff_draws.csv"))
        }
      }
      diff_df <- bind_rows(diff_tbl)
      ev_df   <- imap_dfr(compact(grp_ev), ~ summ(.x) %>% mutate(stage = .y, .before = 1))
      did_df  <- bind_rows(did_tbl)

      ## Write in the brms folder layout ---------------------------------------
      for (src in c("epvs", "vessel", "joint")) {
        cls  <- switch(src, epvs = "EPVS", vessel = "vessel", joint = c("EPVS", "vessel"))
        adir <- sub_dir(reg_dir, "group", src, "all_stages")
        ast  <- file.path(adir, paste(scl, prp, rgn, src, "all_stages", sep = "_"))
        write_csv(fit_tbl %>% filter(class %in% cls), paste0(ast, "_fitted.csv"))
        if (nrow(diff_df) > 0) write_csv(diff_df %>% filter(class %in% cls), paste0(ast, "_model_diff.csv"))
        if (src == "joint") {
          if (nrow(ev_df)  > 0) write_csv(ev_df,  paste0(ast, "_EPVS_minus_vessel.csv"))
          if (nrow(did_df) > 0) write_csv(did_df, paste0(ast, "_EPVS_minus_vessel_diff.csv"))
        }
      }
      message("  subjects: ", paste0(subj_tbl$subjectID, " (", subj_tbl$Stage, ")", collapse = ", "))
      message("  Finished: ", reg_dir)
    }
  }
}
message("\nAll done: ", root_dir)
