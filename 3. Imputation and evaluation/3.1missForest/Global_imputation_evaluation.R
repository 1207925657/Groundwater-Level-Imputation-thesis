if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  if (!nzchar(Sys.setlocale("LC_CTYPE", ".UTF-8"))) stop("A UTF-8 Windows R runtime is required.")
}

read_table <- function(path) as.data.frame(data.table::fread(
  path, encoding = "UTF-8", check.names = FALSE, na.strings = c("", "NA", "NaN")))
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, bom = TRUE, na = "", dateTimeAs = "ISO")
}
parse_flag <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (anyNA(x) || any(!x %in% c("true", "false", "1", "0"))) stop("Invalid removal flags.")
  x %in% c("true", "1")
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

load_well <- function(well, input_dir) {
  data <- read_table(file.path(input_dir, well, "data.csv"))
  audit <- read_table(file.path(input_dir, well, "mask_records.csv"))
  if (length(setdiff(c("date", "target_gwl_original"), names(data))) ||
      length(setdiff(c("date", "target_gwl_original", "target_gwl_observed", "sgi", "was_removed"), names(audit)))) {
    stop("Missing target or mask-record columns.")
  }
  data$date <- as.Date(data$date)
  audit$date <- as.Date(audit$date)
  if (anyNA(data$date) || anyDuplicated(data$date) || !identical(data$date, audit$date)) {
    stop("Model data and audit dates are not aligned.")
  }
  removed <- parse_flag(audit$was_removed)
  expected_removed <- is.finite(audit$target_gwl_observed) & is.finite(audit$sgi) & audit$sgi < -1
  if (!identical(removed, expected_removed)) stop("Removal flags disagree with SGI drought values.")
  expected <- audit$target_gwl_observed
  expected[removed] <- NA_real_
  if (!isTRUE(all.equal(data$target_gwl_original, expected)) ||
      !isTRUE(all.equal(data$target_gwl_original, audit$target_gwl_original))) {
    stop("Masked target differs from the original observations and removal records.")
  }
  list(data = data, audit = audit, removed = removed)
}

model_input <- function(data, features) {
  forbidden <- c("date", "target_gwl_original", "target_gwl_observed", "target_gwl_imputed",
                 "was_removed", "was_imputed", "was_filled", "evaluation_used", "sgi", "is_sgi_drought")
  if (!length(features) || anyDuplicated(features) || any(features %in% forbidden)) stop("Invalid predictor list.")
  if (length(setdiff(features, names(data)))) stop("A configured predictor column is absent.")
  empty <- features[vapply(data[features], function(x) all(is.na(x)), logical(1))]
  if (length(empty)) return(list(all_missing = empty))
  selected <- data[, c("target_gwl_original", features), drop = FALSE]
  if (!all(vapply(selected, is.numeric, logical(1)))) stop("Model inputs must be numeric.")
  if (any(is.infinite(as.matrix(selected)))) stop("Infinite model inputs.")
  constant <- features[vapply(selected[features], function(x) length(unique(na.omit(x))) <= 1L, logical(1))]
  kept <- setdiff(features, constant)
  if (!length(kept)) stop("No varying predictors remain.")
  if (length(unique(na.omit(selected$target_gwl_original))) < 2L) stop("Insufficient observed target variation.")
  list(data = selected[, c("target_gwl_original", kept), drop = FALSE],
       kept = kept, constant = constant, all_missing = character())
}

run_one <- function(index) {
  context <- JOB_CONTEXT
  opts <- context$options
  job <- context$jobs[index, ]
  well <- job$target_well
  checkpoint <- file.path(opts$output_dir, "checkpoints", well, paste0(job$configuration, ".rds"))
  relative <- file.path("series", well, paste0(job$configuration, ".csv"))
  series_file <- file.path(opts$output_dir, relative)
  fingerprint <- paste(context$signature, context$hashes[[well]], job$configuration,
                       job$status, job$feature_columns, sep = "|")
  if (opts$resume && file.exists(checkpoint)) {
    old <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
    if (!is.null(old) && identical(old$fingerprint, fingerprint) && old$metrics$Status == "complete" &&
        file.exists(series_file) && identical(unname(tools::md5sum(series_file)), old$output_md5)) return(old$metrics)
  }
  row <- list(Group = job$group, Configuration = job$configuration, Target_Well = well,
    Status = "failed", Required_Reference_Wells = job$required_reference_wells,
    Available_Reference_Wells = job$available_reference_wells,
    Requested_Features = job$feature_columns, Feature_Names = "", Number_of_Features = NA_integer_,
    Dropped_Constant_Features = "", All_Missing_Features = "",
    Missing_Target_Before = NA_integer_, Missing_Target_After = NA_integer_,
    N_Requested = NA_integer_, N_Evaluated = 0L, Evaluation_Coverage = NA_real_,
    RMSE = NA_real_, MAE = NA_real_, Pearson_r = NA_real_, Pearson_r_status = "not_evaluated",
    Elapsed_Seconds = NA_real_, Output_File = "", Reason = "", Warnings = "")
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  row <- tryCatch(withCallingHandlers({

    if (!identical(WELL_CACHE$well, well)) {
      WELL_CACHE$input <- load_well(well, opts$input_dir)
      WELL_CACHE$well <- well
    }
    input <- WELL_CACHE$input
    row$N_Requested <- sum(input$removed)
    row$Missing_Target_Before <- sum(is.na(input$data$target_gwl_original))
    if (job$status == "insufficient_reference_wells") {
      row$Status <- "skipped_insufficient_references"
      row$Reason <- sprintf("Requires %d references; only %d available.",
                            job$required_reference_wells, job$available_reference_wells)
    } else {
      selected <- model_input(input$data, strsplit(job$feature_columns, ";", fixed = TRUE)[[1]])
      if (length(selected$all_missing)) {
        row$Status <- "skipped_all_missing_predictors"
        row$All_Missing_Features <- paste(selected$all_missing, collapse = ";")
        row$Reason <- "No observations for required predictors; configuration was not reduced."
      } else {
        row$Feature_Names <- paste(selected$kept, collapse = ";")
        row$Number_of_Features <- length(selected$kept)
        row$Dropped_Constant_Features <- paste(selected$constant, collapse = ";")
        RNGkind("Mersenne-Twister", "Inversion", "Rejection")
        set.seed(opts$seed)
        fit <- missForest::missForest(selected$data, ntree = opts$ntree, maxiter = opts$maxiter,
          backend = opts$backend, num.threads = 1L, parallelize = "no", verbose = FALSE)
        predicted <- input$data$target_gwl_original
        missing <- is.na(predicted)
        predicted[missing] <- fit$ximp[, "target_gwl_original"][missing]
        row$Missing_Target_After <- sum(!is.finite(predicted))
        metrics <- calculate_metrics(input$audit$target_gwl_observed[input$removed], predicted[input$removed])
        row[names(metrics)] <- metrics
        row$Evaluation_Coverage <- if (row$N_Requested) row$N_Evaluated / row$N_Requested else NA_real_
        row$Status <- if (!row$N_Requested) "no_evaluation_points" else if (row$N_Evaluated == row$N_Requested) {
          "complete"
        } else "partial"
        output <- data.frame(date = input$data$date, target_gwl_observed = input$audit$target_gwl_observed,
          target_gwl_original = input$data$target_gwl_original, target_gwl_imputed = predicted,
          was_removed = input$removed, was_filled = missing & is.finite(predicted),
          evaluation_used = input$removed & is.finite(input$audit$target_gwl_observed) & is.finite(predicted))
        write_table(output, series_file)
        row$Output_File <- gsub("\\", "/", relative, fixed = TRUE)
      }
    }
    row
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    row$Status <- "failed"
    row$Output_File <- ""
    row$Reason <- conditionMessage(e)
    row
  })
  row$Warnings <- paste(unique(warnings), collapse = " | ")
  row$Elapsed_Seconds <- proc.time()[["elapsed"]] - started
  result <- as.data.frame(row, stringsAsFactors = FALSE)
  dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fingerprint = fingerprint, metrics = result,
    output_md5 = if (nzchar(result$Output_File)) unname(tools::md5sum(series_file)) else NA_character_), checkpoint)
  cat(sprintf("%s %s %s (%.1f s)\n", result$Status, well, job$configuration, row$Elapsed_Seconds))
  flush.console()
  result
}

summarize_metrics <- function(metrics, configs, output_dir) {
  write_table(metrics, file.path(output_dir, "metrics_by_well_configuration.csv"))
  write_table(metrics[startsWith(metrics$Status, "skipped_"), ], file.path(output_dir, "skipped_configurations.csv"))

  completed_sets <- lapply(configs$configuration, function(x) {
    metrics$Target_Well[metrics$Configuration == x & metrics$Status == "complete"]
  })
  common <- sort(Reduce(intersect, completed_sets))
  write_table(data.frame(Target_Well = common), file.path(output_dir, "common_evaluated_wells.csv"))
  average <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  median_value <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  tables <- lapply(c("available_wells", "common_wells"), function(scope) {
    rows <- lapply(seq_len(nrow(configs)), function(i) {
      selected <- metrics[metrics$Configuration == configs$configuration[i], ]
      valid <- selected[selected$Status == "complete", ]
      if (scope == "common_wells") valid <- valid[valid$Target_Well %in% common, ]
      data.frame(Group = configs$group[i], Configuration = configs$configuration[i], Scope = scope,
        Wells_Requested = nrow(selected), Wells_Completed = sum(selected$Status == "complete"),
        Wells_Skipped_References = sum(selected$Status == "skipped_insufficient_references"),
        Wells_Skipped_Predictors = sum(selected$Status == "skipped_all_missing_predictors"),
        Wells_Failed = sum(selected$Status == "failed"), Wells_Partial = sum(selected$Status == "partial"),
        Wells_Without_Evaluation = sum(selected$Status == "no_evaluation_points"),
        Wells_In_Summary = nrow(valid), Wells_With_Defined_Pearson_r = sum(is.finite(valid$Pearson_r)),
        N_Evaluated = sum(valid$N_Evaluated),
        Mean_RMSE = average(valid$RMSE), Median_RMSE = median_value(valid$RMSE),
        Mean_MAE = average(valid$MAE), Median_MAE = median_value(valid$MAE),
        Mean_Pearson_r = average(valid$Pearson_r), Median_Pearson_r = median_value(valid$Pearson_r))
    })
    do.call(rbind, rows)
  })
  write_table(tables[[1]], file.path(output_dir, "metrics_by_configuration.csv"))
  write_table(tables[[2]], file.path(output_dir, "metrics_common_wells.csv"))
}

parse_options <- function(args, here) {
  opts <- list(input_dir = file.path(dirname(dirname(here)), "2.data_preprocessing", "2.2 missForest", "Global"),
    output_dir = file.path(here, "Global"), wells = NULL, configurations = NULL, groups = NULL,
    workers = 4L, ntree = 100L, maxiter = 10L, seed = 42L, backend = "ranger", resume = TRUE, dry_run = FALSE)
  if ("--help" %in% args) {
    cat("Rscript impute_evaluate_global.R [options]\n",
      "--input-dir PATH  --output-dir PATH  --wells Motala_13,Lerum_1\n",
      "--configurations B8_Top08_GWL,B9_Top10_GWL  --groups Baseline,Main_route\n",
      "--workers 4  --ntree 100  --maxiter 10  --seed 42\n",
      "--backend ranger|randomForest  --no-resume  --dry-run\n", sep = "")
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    if (args[i] %in% c("--no-resume", "--dry-run")) {
      if (args[i] == "--no-resume") opts$resume <- FALSE else opts$dry_run <- TRUE
      i <- i + 1L
      next
    }
    name <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
    if (!startsWith(args[i], "--") || !name %in% setdiff(names(opts), c("resume", "dry_run")) ||
        i == length(args)) stop("Invalid option: ", args[i])
    value <- args[i + 1L]
    if (name %in% c("wells", "configurations", "groups")) value <- strsplit(value, ",", fixed = TRUE)[[1]]
    if (name %in% c("workers", "ntree", "maxiter", "seed")) {
      value <- suppressWarnings(as.numeric(value))
      if (!is.finite(value) || value != floor(value) || value > .Machine$integer.max ||
          value < if (name == "seed") 0 else 1) stop("Invalid integer option: ", name)
      value <- as.integer(value)
    }
    opts[[name]] <- value
    i <- i + 2L
  }
  if (!opts$backend %in% c("ranger", "randomForest")) stop("Unsupported backend.")
  opts$input_dir <- normalizePath(opts$input_dir, winslash = "/", mustWork = TRUE)
  opts$output_dir <- normalizePath(opts$output_dir, winslash = "/", mustWork = FALSE)
  opts
}

prepare_jobs <- function(opts) {
  configs <- read_table(file.path(opts$input_dir, "experiment_configurations.csv"))
  manifest <- read_table(file.path(opts$input_dir, "experiment_manifest.csv"))
  if (anyDuplicated(configs$configuration)) stop("Duplicate configuration IDs.")

  for (n in c(8L, 10L)) {
    id <- if (n == 8L) "B8_Top08_GWL" else "B9_Top10_GWL"
    row <- configs[configs$configuration == id, ]
    expected <- paste(sprintf("nearby_%02d_gwl", seq_len(n)), collapse = ";")
    if (nrow(row) != 1L || row$reference_wells != n || row$feature_columns != expected) {
      stop("Regenerate global preprocessing: missing or incorrect configuration ", id)
    }
  }
  jobs <- merge(manifest, configs, by = "configuration", sort = FALSE)
  if (nrow(jobs) != nrow(manifest) || anyDuplicated(jobs[c("target_well", "configuration")])) {
    stop("Configuration/manifest keys are inconsistent.")
  }
  for (pair in list(c("wells", "target_well"), c("configurations", "configuration"), c("groups", "group"))) {
    chosen <- opts[[pair[1]]]
    if (!is.null(chosen)) {
      if (length(setdiff(chosen, jobs[[pair[2]]]))) stop("Unknown selection: ", pair[1])
      jobs <- jobs[jobs[[pair[2]]] %in% chosen, ]
    }
  }
  if (!nrow(jobs)) stop("No matching configurations.")
  ids <- c(jobs$target_well, jobs$configuration)
  if (anyNA(ids) || any(!nzchar(ids)) || any(ids %in% c(".", "..")) || any(basename(ids) != ids) ||
      any(grepl(":", ids, fixed = TRUE)) || any(grepl("\\", ids, fixed = TRUE))) stop("Unsafe output identifier.")
  if (any(!jobs$status %in% c("ready", "insufficient_reference_wells")) ||
      any(jobs$required_reference_wells != jobs$reference_wells) ||
      any((jobs$status == "ready") != (jobs$available_reference_wells >= jobs$reference_wells))) {
    stop("Inconsistent preprocessing availability records.")
  }
  jobs$planned_status <- ifelse(jobs$status != "ready", "skipped_insufficient_references",
    ifelse(!is.na(jobs$all_missing_predictor_columns) & nzchar(jobs$all_missing_predictor_columns),
           "skipped_all_missing_predictors", "pending"))
  jobs <- jobs[order(jobs$target_well, -jobs$number_of_features, jobs$configuration), ]
  list(jobs = jobs, configs = configs[configs$configuration %in% jobs$configuration, ])
}

main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), winslash = "/")
  opts <- parse_options(commandArgs(trailingOnly = TRUE), dirname(script))
  if (is.null(opts)) return(invisible(NULL))
  for (p in c("missForest", "data.table")) if (!requireNamespace(p, quietly = TRUE)) stop("Install required R package: ", p)
  prepared <- prepare_jobs(opts)
  jobs <- prepared$jobs
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_table(jobs, file.path(opts$output_dir, if (opts$dry_run) "planned_jobs.csv" else "run_manifest.csv"))
  cat(sprintf("Selected %d pairs (%d wells, %d configurations).\n", nrow(jobs),
              length(unique(jobs$target_well)), nrow(prepared$configs)))
  print(table(jobs$planned_status))
  if (opts$dry_run) return(invisible(jobs))
  wells <- unique(jobs$target_well)
  hashes <- setNames(lapply(wells, function(well) paste(unname(tools::md5sum(
    file.path(opts$input_dir, well, c("data.csv", "mask_records.csv")))), collapse = ":")), wells)
  if (any(grepl("NA", unlist(hashes), fixed = TRUE))) stop("Missing well input files.")
  versions <- vapply(c("missForest", "ranger", "randomForest", "data.table"),
                     function(p) as.character(utils::packageVersion(p)), character(1))
  signature <- paste("global-v1", paste(unname(tools::md5sum(c(script,
    file.path(opts$input_dir, c("experiment_configurations.csv", "experiment_manifest.csv"))))), collapse = ":"),
    R.version.string, paste(versions, collapse = "/"), opts$ntree, opts$maxiter, opts$seed, opts$backend, sep = "|")
  JOB_CONTEXT <<- list(options = opts, jobs = jobs, hashes = hashes, signature = signature)
  WELL_CACHE <<- new.env(parent = emptyenv())
  writeLines(c(paste(names(opts), vapply(opts, function(x) paste(x, collapse = ","), character(1)), sep = " = "),
    paste(names(versions), versions, sep = " = "), paste("signature =", signature)),
    file.path(opts$output_dir, "run_settings.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(opts$output_dir, "sessionInfo.txt"))
  if (opts$workers == 1L) {
    results <- lapply(seq_len(nrow(jobs)), run_one)
  } else {
    cluster <- parallel::makePSOCKcluster(min(opts$workers, nrow(jobs)),
                                         outfile = file.path(opts$output_dir, "workers.log"))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(path) source(path, encoding = "UTF-8"), script)
    parallel::clusterExport(cluster, "JOB_CONTEXT", envir = .GlobalEnv)
    parallel::clusterEvalQ(cluster, { data.table::setDTthreads(1L); WELL_CACHE <- new.env(parent = emptyenv()); NULL })
    results <- parallel::parLapplyLB(cluster, seq_len(nrow(jobs)), run_one)
  }
  metrics <- do.call(rbind, results)
  metrics <- metrics[order(match(metrics$Configuration, prepared$configs$configuration), metrics$Target_Well), ]
  summarize_metrics(metrics, prepared$configs, opts$output_dir)
  cat(sprintf("Finished: %d complete, %d skipped, %d failed/partial. Results: %s\n",
    sum(metrics$Status == "complete"), sum(startsWith(metrics$Status, "skipped_")),
    sum(metrics$Status %in% c("failed", "partial")), opts$output_dir))
  if (any(metrics$Status %in% c("failed", "partial"))) stop("Inspect failed/partial cases in the saved metrics and resume.")
  invisible(metrics)
}

if (sys.nframe() == 0L) main()
