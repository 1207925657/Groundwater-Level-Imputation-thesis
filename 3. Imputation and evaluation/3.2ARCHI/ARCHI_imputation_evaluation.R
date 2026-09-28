if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  if (!nzchar(Sys.setlocale("LC_CTYPE", ".UTF-8"))) stop("Use a UTF-8 Windows R runtime (R >= 4.2).")
}
read_table <- function(path) as.data.frame(data.table::fread(
  path, encoding = "UTF-8", check.names = FALSE, na.strings = c("", "NA", "NaN")))
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, bom = TRUE, na = "", dateTimeAs = "ISO")
}
csv_date <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  unique_dates <- unique(as.character(x))
  dates <- as.Date(unique_dates, format = "%m/%d/%Y")
  if (anyNA(dates)) stop("Invalid preprocessing date; expected MM/DD/YYYY.")
  dates[match(x, unique_dates)]
}
valid_id <- function(x) {
  !is.na(x) & nzchar(x) & !x %in% c(".", "..") & basename(x) == x &
    !grepl("\\", x, fixed = TRUE) & !grepl(":", x, fixed = TRUE)
}
calculate_metrics <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  n <- length(observed)
  reason <- if (n < 2L) "fewer_than_two_pairs" else if (sd(observed) == 0) {
    "constant_observed"
  } else if (sd(predicted) == 0) "constant_predicted" else "defined"
  list(N_Evaluated = n,
    RMSE = if (n) sqrt(mean((predicted - observed)^2)) else NA_real_,
    MAE = if (n) mean(abs(predicted - observed)) else NA_real_,
    Pearson_r = if (reason == "defined") cor(observed, predicted, method = "pearson") else NA_real_,
    Pearson_r_status = reason)
}

prepare_target <- function(gwl, job, masks) {
  target <- gwl[gwl$site_no == job$target_well, c("date", "value")]
  target <- target[order(target$date), ]
  expected_dates <- seq.Date(job$target_first_date, job$eval_end, by = "day")
  if (!identical(target$date, expected_dates)) stop("Target calendar disagrees with the manifest.")
  removed <- masks[masks$site_no == job$target_well, ]
  position <- match(removed$date, target$date)
  if (nrow(removed) != job$masked_target_rows || anyNA(position) || anyDuplicated(position) ||
      any(!is.na(target$value[position]))) stop("Target values and removal records disagree.")
  truth <- target$value
  truth[position] <- removed$original_value
  flags <- rep(FALSE, nrow(target))
  flags[position] <- TRUE
  data.frame(date = target$date, target_gwl_observed = truth,
    target_gwl_original = target$value, target_gwl_imputed = target$value,
    was_removed = flags, was_filled = FALSE, evaluation_used = FALSE,
    archi_time_range = NA_character_)
}

run_target <- function(job, context) {
  opts <- context$options
  well <- job$target_well
  source_file <- file.path(opts$input_dir, "Scenario_A", job$output_file)
  checkpoint <- file.path(opts$output_dir, "checkpoints", paste0(well, ".rds"))
  source_hash <- unname(tools::md5sum(source_file))
  fingerprint <- paste(context$signature, source_hash, job$random_seed, sep = "|")
  if (opts$resume && !is.na(source_hash) && file.exists(checkpoint)) {
    old <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
    if (!is.null(old) && identical(old$fingerprint, fingerprint) && old$metrics$Status != "failed") {
      files <- file.path(opts$output_dir, names(old$output_md5))
      if (length(files) && all(file.exists(files)) &&
          identical(unname(tools::md5sum(files)), unname(old$output_md5))) return(old$metrics)
    }
  }
  row <- list(Target_Well = well, Status = "failed", Target_Start = as.character(job$target_first_date),
    Target_End = as.character(job$eval_end), Random_Seed = job$random_seed,
    Maximum_Reference_Wells = opts$n_refwl, Reference_Wells_Used = 0L,
    Dropped_Short_Reference_Wells = NA_integer_, Candidate_Wells_In_Grid = NA_integer_,
    Grid_Rows_Before_Trim = NA_integer_, Grid_Rows_After_Trim = NA_integer_,
    Missing_Target_Before = NA_integer_, Missing_Target_After = NA_integer_,
    N_Requested = job$masked_target_rows, N_Evaluated = 0L, Evaluation_Coverage = 0,
    RMSE = NA_real_, MAE = NA_real_, Pearson_r = NA_real_, Pearson_r_status = "not_evaluated",
    Regression_Model = "", Model_Notes = "", Reason = "", Warnings = "",
    Elapsed_Seconds = NA_real_, Output_File = "")
  outputs <- character()
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  cat("Running ARCHI Scenario A:", well, "\n")
  row <- tryCatch(withCallingHandlers({
    if (is.na(source_hash)) stop("Missing Scenario A input CSV.")
    gwl <- read_table(source_file)
    if (!identical(names(gwl), c("site_no", "date", "value")) || !is.numeric(gwl$value)) {
      stop("ARCHI input must have exactly site_no,date,value, with numeric values.")
    }
    if (anyNA(gwl$site_no) || any(is.infinite(gwl$value)) || nrow(gwl) != job$total_rows ||
        !setequal(unique(gwl$site_no), context$metadata$site_no)) stop("Invalid full-network input.")
    gwl$date <- csv_date(gwl$date)
    series <- prepare_target(gwl, job, context$masks)
    row$Missing_Target_Before <- sum(is.na(series$target_gwl_original))
    end_dates <- context$ends$eval_end[match(gwl$site_no, context$ends$site_no)]
    if (anyNA(end_dates)) stop("A reference well has no original record end date.")
    eligible <- context$ends$site_no[context$ends$eval_end >= job$eval_end]
    eligible <- unique(c(well, eligible))
    row$Dropped_Short_Reference_Wells <- length(setdiff(unique(gwl$site_no), eligible))
    gwl <- gwl[gwl$site_no %in% eligible & gwl$date <= end_dates &
                 gwl$date >= job$target_first_date & gwl$date <= job$eval_end, ]
    if (anyDuplicated(gwl[c("site_no", "date")])) stop("Duplicate site/date input rows.")

    wide <- data.table::dcast(data.table::as.data.table(gwl), date ~ site_no, value.var = "value")
    data.table::setnames(wide, "date", "timestep")
    data.table::setcolorder(wide, c("timestep", sort(unique(gwl$site_no))))
    grid <- as.data.frame(wide)
    grid$timestep <- as.Date(grid$timestep)
    if (!identical(grid$timestep, series$date)) stop("Daily grid lost target dates.")
    rm(gwl, wide, end_dates)
    row$Grid_Rows_Before_Trim <- nrow(grid)
    grid <- ARCHI::trim_grid(grid, data_thresh = 0.35, time_thresh = 0, rm_nzv = TRUE)
    row$Grid_Rows_After_Trim <- nrow(grid)
    row$Candidate_Wells_In_Grid <- sum(!names(grid) %in% c("timestep", well))
    stats <- data.frame(site_no = character(), n_refs = integer(), mod_notes = character())
    refs <- character()
    if (!well %in% names(grid)) {
      row$Status <- "target_excluded"
      row$Reason <- "Target excluded by completeness or near-zero-variance trimming."
    } else if (ncol(grid) < 3L) {
      row$Status <- "not_imputed"
      row$Reason <- "No reference site remains in the grid."
    } else {
      metadata <- context$metadata[context$metadata$site_no %in% names(grid), ]
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      set.seed(job$random_seed)
      fit <- ARCHI::impute_grid(grid, model = "ridge", error_method = "NSE", error_thresh = 0,
        n_refwl = opts$n_refwl, p_per_n = NA_real_, r_cutoff = NA_real_,
        add_std_means = FALSE, sites = metadata, sites_crs = "WGS84", d_cutoff = 50,
        group_sites = FALSE, relax = 0.1, final_pass = TRUE, rnd = 2,
        bootstrap_PI = FALSE, cv_lambda = "lambda.min", nfolds = 10, verbose = FALSE)
      stats <- as.data.frame(fit$model_stats)
      stats <- stats[stats$site_no == well, , drop = FALSE]
      if (nrow(stats) != 1L) stop("Missing or duplicate target model statistics.")
      row$Model_Notes <- if (is.na(stats$mod_notes[1])) "" else as.character(stats$mod_notes[1])
      row$Regression_Model <- if (is.na(stats$regression_model[1])) "" else as.character(stats$regression_model[1])
      refs <- as.character(fit$refs[[well]])
      refs <- refs[!is.na(refs) & nzchar(refs)]
      if (anyDuplicated(refs) || any(!refs %in% metadata$site_no) || well %in% refs) stop("Invalid selected reference list.")
      row$Reference_Wells_Used <- if (is.na(stats$n_refs[1])) 0L else as.integer(stats$n_refs[1])
      if (row$Reference_Wells_Used != length(refs)) stop("Reference count and names disagree.")
      filled_grid <- fit$imputed_grid
      if (!well %in% names(filled_grid) || anyDuplicated(filled_grid$timestep)) stop("Invalid imputed grid.")
      predictions <- filled_grid[[well]][match(series$date, as.Date(filled_grid$timestep))]
      missing <- is.na(series$target_gwl_original)
      series$target_gwl_imputed[missing] <- predictions[missing]
      series$was_filled <- missing & is.finite(series$target_gwl_imputed)
      long <- as.data.frame(fit$out_long)
      long <- long[!is.na(long$site_no) & long$site_no == well & long$value_type %in% "Imputed", ]
      long$timestep <- as.Date(long$timestep)
      if (anyDuplicated(long$timestep)) stop("Duplicate ARCHI imputation dates.")
      positions <- match(series$date[series$was_filled], long$timestep)

      if (length(positions) && (anyNA(positions) || !isTRUE(all.equal(series$target_gwl_imputed[series$was_filled],
                                               long$value[positions], check.attributes = FALSE)))) {
        stop("Filled values disagree with ARCHI's Imputed records.")
      }
      series$archi_time_range[series$was_filled] <- long$time_range[positions]
      row$Status <- "evaluated"
    }
    series$evaluation_used <- series$was_removed & series$was_filled &
      is.finite(series$target_gwl_observed) & is.finite(series$target_gwl_imputed)
    metrics <- calculate_metrics(series$target_gwl_observed[series$evaluation_used],
                                 series$target_gwl_imputed[series$evaluation_used])
    row[names(metrics)] <- metrics
    row$Missing_Target_After <- sum(!is.finite(series$target_gwl_imputed))
    row$Evaluation_Coverage <- if (row$N_Requested) row$N_Evaluated / row$N_Requested else NA_real_
    if (row$Status == "evaluated") {
      row$Status <- if (!row$N_Requested) "no_evaluation_points" else if (!row$N_Evaluated) {
        "not_imputed"
      } else if (row$N_Evaluated == row$N_Requested) "complete" else "partial"
      if (row$Status == "not_imputed") row$Reason <- row$Model_Notes
    }
    ref_table <- data.frame(Target_Well = rep(well, length(refs)),
                           Reference_Rank = seq_along(refs), Reference_Well = refs)
    outputs <- file.path(c("series", "model_stats", "reference_wells"), paste0(well, ".csv"))
    write_table(series, file.path(opts$output_dir, outputs[1]))
    write_table(stats, file.path(opts$output_dir, outputs[2]))
    write_table(ref_table, file.path(opts$output_dir, outputs[3]))
    row$Output_File <- gsub("\\", "/", outputs[1], fixed = TRUE)
    row
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    row$Status <- "failed"
    row$Reason <- conditionMessage(e)
    row$Output_File <- ""
    row
  })
  row$Warnings <- paste(unique(warnings), collapse = " | ")
  row$Elapsed_Seconds <- proc.time()[["elapsed"]] - started
  result <- as.data.frame(row, stringsAsFactors = FALSE)
  dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
  hashes <- if (row$Status == "failed") character() else setNames(
    unname(tools::md5sum(file.path(opts$output_dir, outputs))), outputs)
  saveRDS(list(fingerprint = fingerprint, metrics = result, output_md5 = hashes), checkpoint)
  cat(sprintf("%s %s: %d/%d removed values evaluated; %d references (%.1f s).\n", well,
              row$Status, row$N_Evaluated, row$N_Requested, row$Reference_Wells_Used, row$Elapsed_Seconds))
  result
}

save_summary <- function(metrics, output_dir) {
  write_table(metrics, file.path(output_dir, "metrics_by_well.csv"))
  mean_finite <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  median_finite <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  rows <- lapply(c("complete_wells", "wells_with_any_predictions"), function(scope) {
    valid <- metrics[metrics$Status %in% if (scope == "complete_wells") "complete" else c("complete", "partial"), ]
    data.frame(Scope = scope, Wells_Requested = nrow(metrics), Wells_In_Summary = nrow(valid),
      Wells_Complete = sum(metrics$Status == "complete"), Wells_Partial = sum(metrics$Status == "partial"),
      Wells_Not_Imputed = sum(metrics$Status == "not_imputed"),
      Wells_Excluded = sum(metrics$Status == "target_excluded"), Wells_Failed = sum(metrics$Status == "failed"),
      Wells_Without_Evaluation = sum(metrics$Status == "no_evaluation_points"),
      Wells_With_Defined_Pearson_r = sum(is.finite(valid$Pearson_r)),
      N_Requested_In_Summary = sum(valid$N_Requested), N_Evaluated = sum(valid$N_Evaluated),
      Mean_RMSE = mean_finite(valid$RMSE), Median_RMSE = median_finite(valid$RMSE),
      Mean_MAE = mean_finite(valid$MAE), Median_MAE = median_finite(valid$MAE),
      Mean_Pearson_r = mean_finite(valid$Pearson_r), Median_Pearson_r = median_finite(valid$Pearson_r))
  })
  write_table(do.call(rbind, rows), file.path(output_dir, "metrics_summary.csv"))
}

parse_options <- function(args, here) {
  opts <- list(input_dir = file.path(dirname(dirname(here)), "2.data_preprocessing", "2.3 ARCHI"),
               output_dir = NULL, wells = NULL, n_refwl = 10L, seed = 123L, resume = TRUE, dry_run = FALSE)
  if ("--help" %in% args) {
    cat("Rscript impute_evaluate_archi.R [--wells Motala_13,Lerum_1] [--n-refwl 10]\n",
        "  [--input-dir PATH] [--output-dir PATH] [--seed 123] [--no-resume] [--dry-run]\n")
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    if (args[i] %in% c("--no-resume", "--dry-run")) {
      if (args[i] == "--no-resume") opts$resume <- FALSE else opts$dry_run <- TRUE
      i <- i + 1L
      next
    }
    key <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
    if (!startsWith(args[i], "--") || !key %in% setdiff(names(opts), c("resume", "dry_run")) ||
        i == length(args)) stop("Invalid option: ", args[i])
    value <- args[i + 1L]
    if (key == "wells") value <- strsplit(value, ",", fixed = TRUE)[[1]]
    if (key %in% c("n_refwl", "seed")) {
      value <- suppressWarnings(as.numeric(value))
      if (!is.finite(value) || value != floor(value) || value > .Machine$integer.max - 100000 ||
          value < if (key == "seed") 0 else 1) stop("Invalid integer option: ", key)
      value <- as.integer(value)
    }
    opts[[key]] <- value
    i <- i + 2L
  }
  opts$input_dir <- normalizePath(opts$input_dir, winslash = "/", mustWork = TRUE)
  if (is.null(opts$output_dir)) opts$output_dir <- file.path(here, "Scenario_A", paste0("nref_", opts$n_refwl))
  opts$output_dir <- normalizePath(opts$output_dir, winslash = "/", mustWork = FALSE)
  opts
}

main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), winslash = "/")
  opts <- parse_options(commandArgs(trailingOnly = TRUE), dirname(script))
  if (is.null(opts)) return(invisible(NULL))
  for (p in c("ARCHI", "data.table")) if (!requireNamespace(p, quietly = TRUE)) stop("Install required R package: ", p)
  files <- file.path(opts$input_dir, c("Scenario_A/ARCHI_1masked_summary.csv", "Scenario_A/site_metadata.csv",
                                      "well_summary.csv", "mask_records.csv"))
  jobs <- read_table(files[1])
  if (anyDuplicated(jobs$target_well) || any(!valid_id(jobs$target_well)) || any(!valid_id(jobs$output_file))) {
    stop("Invalid Scenario A manifest identifiers.")
  }
  jobs$random_seed <- opts$seed + seq_len(nrow(jobs))
  jobs$target_first_date <- csv_date(jobs$target_first_date)
  jobs$eval_end <- csv_date(jobs$eval_end)
  if (!is.null(opts$wells)) {
    if (length(setdiff(opts$wells, jobs$target_well))) stop("Unknown target well selection.")
    jobs <- jobs[jobs$target_well %in% opts$wells, ]
  }
  if (!nrow(jobs)) stop("No target wells selected.")
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_table(jobs, file.path(opts$output_dir, if (opts$dry_run) "planned_targets.csv" else "run_manifest.csv"))
  if (opts$dry_run) {
    cat(nrow(jobs), "target wells planned; no ARCHI model was fitted.\n")
    return(invisible(jobs))
  }
  metadata <- read_table(files[2])
  ends <- read_table(files[3])
  masks <- read_table(files[4])
  if (!identical(names(metadata), c("site_no", "latitude", "longitude", "group")) ||
      anyDuplicated(metadata$site_no) || anyDuplicated(ends$site_no) ||
      !setequal(metadata$site_no, ends$site_no) || any(!jobs$target_well %in% metadata$site_no) ||
      any(!is.finite(metadata$latitude)) || any(!is.finite(metadata$longitude))) stop("Invalid site metadata.")
  ends$eval_end <- csv_date(ends$eval_end)
  if (!identical(jobs$eval_end, ends$eval_end[match(jobs$target_well, ends$site_no)])) stop("Conflicting record end dates.")
  masks$date <- csv_date(masks$date)
  flags <- tolower(as.character(masks$was_removed))
  if (anyNA(flags) || any(!flags %in% c("true", "1")) || any(!is.finite(masks$original_value)) ||
      any(!is.finite(masks$sgi) | masks$sgi >= -1) || anyDuplicated(masks[c("site_no", "date")])) {
    stop("Invalid artificial-removal records.")
  }
  versions <- vapply(c("ARCHI", "glmnet", "data.table", "caret", "terra", "dplyr"),
    function(p) as.character(utils::packageVersion(p)), character(1))
  remote_sha <- utils::packageDescription("ARCHI")$RemoteSha
  signature <- paste("archi-v1", paste(unname(tools::md5sum(c(script, files))), collapse = ":"),
    R.version.string, paste(versions, collapse = "/"), paste(remote_sha, collapse = ""), opts$n_refwl, sep = "|")
  context <- list(options = opts, metadata = metadata, ends = ends, masks = masks, signature = signature)
  writeLines(c(paste(names(opts), vapply(opts, function(x) paste(x, collapse = ","), character(1)), sep = " = "),
    paste(names(versions), versions, sep = " = "), paste("ARCHI_RemoteSha =", remote_sha),
    "daily; ridge; NSE acceptance threshold 0; data_thresh 0.35; distance 50 km; add_std_means FALSE",
    "p_per_n NA; r_cutoff NA; group_sites FALSE; drop references ending before target TRUE",
    "relax 0.1; final_pass TRUE; rnd 2; bootstrap_PI FALSE; cv_lambda lambda.min; nfolds 10",
    paste("signature =", signature)), file.path(opts$output_dir, "run_settings.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(opts$output_dir, "sessionInfo.txt"))
  results <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    results[[i]] <- run_target(jobs[i, ], context)
    save_summary(do.call(rbind, results[seq_len(i)]), opts$output_dir)
    invisible(gc())
  }
  metrics <- do.call(rbind, results)
  cat("Saved", nrow(metrics), "target records to", opts$output_dir, "\n")
  if (any(metrics$Status == "failed")) stop("Some cases failed; inspect metrics_by_well.csv and resume.")
  invisible(metrics)
}

if (sys.nframe() == 0L) main()
