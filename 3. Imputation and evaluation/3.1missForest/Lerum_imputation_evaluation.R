read_table <- function(path) {
  as.data.frame(data.table::fread(path, encoding = "UTF-8", check.names = FALSE,
                                 na.strings = c("", "NA", "NaN")))
}

write_table <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data, path, bom = TRUE, na = "", dateTimeAs = "ISO")
}

parse_flag <- function(x) {
  text <- tolower(trimws(as.character(x)))
  if (anyNA(text) || any(!text %in% c("true", "false", "1", "0"))) {
    stop("Invalid or missing removal flag.")
  }
  text %in% c("true", "1")
}

calculate_metrics <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  n <- length(observed)
  pearson_status <- if (n < 2L) "fewer_than_two_pairs" else if (sd(observed) == 0) {
    "constant_observed"
  } else if (sd(predicted) == 0) "constant_predicted" else "defined"
  list(N_Evaluated = n,
       RMSE = if (n) sqrt(mean((predicted - observed)^2)) else NA_real_,
       MAE = if (n) mean(abs(predicted - observed)) else NA_real_,
       Pearson_r = if (pearson_status == "defined") {
         cor(observed, predicted, method = "pearson")
       } else NA_real_,
       Pearson_r_status = pearson_status)
}

load_well <- function(well, input_dir) {
  data <- read_table(file.path(input_dir, "data", paste0(well, ".csv")))
  audit <- read_table(file.path(input_dir, "mask_records", paste0(well, ".csv")))
  data$date <- as.Date(data$date)
  audit$date <- as.Date(audit$date)
  stopifnot(!anyNA(data$date), !anyDuplicated(data$date),
            identical(data$date, audit$date),
            isTRUE(all.equal(data$target_gwl_original, audit$target_gwl_original)))
  removed <- parse_flag(audit$was_removed)
  if (!identical(removed, is.finite(audit$target_gwl_observed) & audit$sgi < -1)) {
    stop("Removal flags disagree with SGI drought values: ", well)
  }
  if (any(removed & !is.na(data$target_gwl_original))) stop("Unmasked held-out target values.")
  observed <- !is.na(data$target_gwl_original)
  if (!isTRUE(all.equal(data$target_gwl_original[observed], audit$target_gwl_observed[observed]))) {
    stop("Observed target values changed before imputation: ", well)
  }
  list(data = data, audit = audit, removed = removed)
}

model_input <- function(data, features) {
  forbidden <- c("date", "target_gwl_original", "target_gwl_observed",
                 "target_gwl_imputed", "was_removed", "was_imputed", "sgi", "is_sgi_drought")
  if (!length(features) || anyDuplicated(features) || any(features %in% forbidden)) {
    stop("Invalid configuration or a forbidden predictor column.")
  }
  if (length(setdiff(features, names(data)))) stop("A configured predictor column is absent.")
  selected <- data[, c("target_gwl_original", features), drop = FALSE]
  if (!all(vapply(selected, is.numeric, logical(1)))) stop("Every model column must be numeric.")
  if (any(vapply(selected, function(x) any(is.infinite(x)), logical(1)))) stop("Infinite model input.")
  drop <- features[vapply(selected[features], function(x) length(unique(x[!is.na(x)])) <= 1L, logical(1))]
  kept <- setdiff(features, drop)
  if (!length(kept)) stop("No varying, observed predictors remain.")
  if (length(unique(na.omit(selected$target_gwl_original))) < 2L) stop("Insufficient observed target variation.")
  list(data = selected[, c("target_gwl_original", kept), drop = FALSE], kept = kept, dropped = drop)
}

run_one <- function(index) {
  context <- JOB_CONTEXT
  opts <- context$options
  job <- context$jobs[index, ]
  key <- paste(job$target_well, job$configuration, sep = "__")
  checkpoint <- file.path(opts$output_dir, "checkpoints", paste0(key, ".rds"))
  series_relative <- file.path("series", job$experiment, paste0(key, ".csv"))
  series_file <- file.path(opts$output_dir, series_relative)
  fingerprint <- paste(context$signature, context$hashes[[job$target_well]],
                       job$experiment, job$configuration, job$feature_columns, sep = "|")
  if (opts$resume && file.exists(checkpoint)) {
    old <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
    if (!is.null(old) && identical(old$fingerprint, fingerprint) &&
        identical(old$metrics$Status, "complete") && file.exists(series_file) &&
        identical(unname(tools::md5sum(series_file)), old$output_md5)) return(old$metrics)
  }
  row <- list(Experiment = job$experiment, Configuration = job$configuration,
              Target_Well = job$target_well, Status = "failed",
              Requested_Features = job$feature_columns, Feature_Names = "",
              Dropped_Features = "", Number_of_Features = NA_integer_,
              Missing_Target_Before = NA_integer_, Missing_Target_After = NA_integer_,
              N_Requested = NA_integer_, N_Evaluated = 0L, Evaluation_Coverage = NA_real_,
              RMSE = NA_real_, MAE = NA_real_, Pearson_r = NA_real_, Pearson_r_status = "not_evaluated",
              Elapsed_Seconds = NA_real_, Output_File = gsub("\\\\", "/", series_relative),
              Error = "", Warnings = "")
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  row <- tryCatch(withCallingHandlers({
    well <- job$target_well
    if (!exists(well, WELL_CACHE, inherits = FALSE)) {
      assign(well, load_well(well, opts$input_dir), envir = WELL_CACHE)
    }
    input <- get(well, WELL_CACHE, inherits = FALSE)
    features <- strsplit(job$feature_columns, ";", fixed = TRUE)[[1]]
    selected <- model_input(input$data, features)
    row$Feature_Names <- paste(selected$kept, collapse = ";")
    row$Dropped_Features <- paste(selected$dropped, collapse = ";")
    row$Number_of_Features <- length(selected$kept)
    row$Missing_Target_Before <- sum(is.na(input$data$target_gwl_original))
    row$N_Requested <- sum(input$removed)

    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(opts$seed)
    fit <- missForest::missForest(selected$data, maxiter = opts$maxiter, ntree = opts$ntree,
                                   backend = opts$backend, num.threads = 1L,
                                   parallelize = "no", verbose = FALSE)
    predicted <- input$data$target_gwl_original
    missing <- is.na(predicted)
    predicted[missing] <- fit$ximp[, "target_gwl_original"][missing]
    if (length(predicted) != nrow(input$data)) stop("Unexpected imputed target length.")
    row$Missing_Target_After <- sum(!is.finite(predicted))
    metrics <- calculate_metrics(input$audit$target_gwl_observed[input$removed], predicted[input$removed])
    row[names(metrics)] <- metrics
    row$Evaluation_Coverage <- if (row$N_Requested) row$N_Evaluated / row$N_Requested else NA_real_
    row$Status <- if (!row$N_Requested) "no_evaluation_points" else if (row$N_Evaluated == row$N_Requested) {
      "complete"
    } else "partial"
    output <- data.frame(date = input$data$date,
                         target_gwl_observed = input$audit$target_gwl_observed,
                         target_gwl_original = input$data$target_gwl_original,
                         target_gwl_imputed = predicted,
                         was_removed = input$removed,
                         was_filled = missing & is.finite(predicted),
                         evaluation_used = input$removed & is.finite(input$audit$target_gwl_observed) & is.finite(predicted))
    write_table(output, series_file)
    row
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    row$Status <- "failed"
    row$Error <- conditionMessage(e)
    row
  })
  row$Elapsed_Seconds <- proc.time()[["elapsed"]] - started
  row$Warnings <- paste(unique(warnings), collapse = " | ")
  result <- as.data.frame(row, stringsAsFactors = FALSE)
  dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fingerprint = fingerprint, metrics = result,
               output_md5 = if (file.exists(series_file)) unname(tools::md5sum(series_file)) else NA_character_), checkpoint)
  cat(sprintf("%s %s %s (%.1f s)\n", result$Status, job$target_well,
              job$configuration, result$Elapsed_Seconds))
  flush.console()
  result
}

summarize_metrics <- function(metrics, jobs, output_dir) {
  write_table(metrics, file.path(output_dir, "metrics_by_well_configuration.csv"))
  configurations <- unique(jobs$configuration)
  average <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  summaries <- lapply(configurations, function(name) {
    rows <- metrics[metrics$Configuration == name, ]
    valid <- rows[rows$Status == "complete", ]
    data.frame(Experiment = rows$Experiment[1], Configuration = name,
               Wells_Requested = nrow(rows), Wells_Completed = nrow(valid),
               Wells_With_Defined_Pearson_r = sum(is.finite(valid$Pearson_r)),
               Mean_RMSE = average(valid$RMSE), Mean_MAE = average(valid$MAE),
               Mean_Pearson_r = average(valid$Pearson_r))
  })
  write_table(do.call(rbind, summaries), file.path(output_dir, "metrics_by_configuration.csv"))
}

parse_options <- function(args, here) {
  repository <- dirname(dirname(here))
  opts <- list(input_dir = file.path(repository, "2.data_preprocessing", "2.2 missForest", "Lerum"),
               output_dir = file.path(here, "Lerum"), wells = NULL, configurations = NULL,
               experiments = NULL, workers = 4L, ntree = 100L, maxiter = 10L,
               seed = 42L, backend = "ranger", resume = TRUE)
  if ("--help" %in% args) {
    cat("Rscript impute_evaluate_lerum.R [options]\n",
        "--input-dir PATH  --output-dir PATH\n",
        "--wells Lerum_1,Lerum_3  --configurations T1_day_of_year,T2_month\n",
        "--experiments Experiment_1_Time,Experiment_2_NP\n",
        "--workers 4  --ntree 100  --maxiter 10  --seed 42\n",
        "--backend ranger|randomForest  --no-resume\n", sep = "")
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    name <- sub("^--", "", args[i])
    if (name == "no-resume") {
      opts$resume <- FALSE
      i <- i + 1L
      next
    }
    name <- gsub("-", "_", name, fixed = TRUE)
    if (!name %in% names(opts) || name == "resume" || i == length(args)) stop("Invalid option: ", args[i])
    value <- args[i + 1L]
    if (name %in% c("wells", "configurations", "experiments")) value <- strsplit(value, ",", fixed = TRUE)[[1]]
    if (name %in% c("workers", "ntree", "maxiter", "seed")) {
      number <- suppressWarnings(as.numeric(value))
      if (!is.finite(number) || number != floor(number) || number < if (name == "seed") 0 else 1) {
        stop("Invalid integer option: ", name)
      }
      value <- as.integer(number)
    }
    opts[[name]] <- value
    i <- i + 2L
  }
  if (!opts$backend %in% c("ranger", "randomForest")) stop("Unsupported forest backend.")
  opts$input_dir <- normalizePath(opts$input_dir, winslash = "/", mustWork = TRUE)
  opts$output_dir <- normalizePath(opts$output_dir, winslash = "/", mustWork = FALSE)
  opts
}

main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), winslash = "/")
  opts <- parse_options(commandArgs(trailingOnly = TRUE), dirname(script))
  if (is.null(opts)) return(invisible(NULL))
  for (package in c("missForest", "data.table")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install required R package: ", package)
  }
  config_path <- file.path(opts$input_dir, "experiment_configurations.csv")
  manifest_path <- file.path(opts$input_dir, "experiment_manifest.csv")
  configs <- read_table(config_path)
  manifest <- read_table(manifest_path)
  if (anyDuplicated(configs$configuration)) stop("Duplicate configuration IDs.")
  jobs <- merge(manifest, configs, by = c("experiment", "configuration"), sort = FALSE)
  if (nrow(jobs) != nrow(manifest) || anyDuplicated(jobs[c("target_well", "configuration")])) {
    stop("Configuration/manifest keys are inconsistent.")
  }
  for (pair in list(c("wells", "target_well"), c("configurations", "configuration"), c("experiments", "experiment"))) {
    selected <- opts[[pair[1]]]
    if (!is.null(selected)) {
      if (length(setdiff(selected, jobs[[pair[2]]]))) stop("Unknown selection: ", pair[1])
      jobs <- jobs[jobs[[pair[2]]] %in% selected, ]
    }
  }
  if (!nrow(jobs)) stop("No matching experiments.")
  if (any(!grepl("^[A-Za-z0-9_]+$", jobs$target_well)) ||
      any(!grepl("^[A-Za-z0-9_]+$", jobs$configuration))) stop("Unsafe output identifier.")
  if (any(jobs$status != "ready")) stop("Some selected preprocessing configurations are unavailable.")

  jobs <- jobs[order(-jobs$number_of_features, jobs$target_well, jobs$configuration), ]
  wells <- sort(unique(jobs$target_well))
  hashes <- setNames(lapply(wells, function(well) {
    paste(unname(tools::md5sum(c(file.path(opts$input_dir, "data", paste0(well, ".csv")),
                                file.path(opts$input_dir, "mask_records", paste0(well, ".csv"))))), collapse = ":")
  }), wells)
  if (any(grepl("NA", unlist(hashes), fixed = TRUE))) stop("Missing well input file.")
  versions <- vapply(c("missForest", "ranger", "randomForest", "data.table"),
                     function(package) as.character(utils::packageVersion(package)), character(1))
  signature <- paste("lerum-v1", unname(tools::md5sum(script)), R.version.string,
                      paste(versions, collapse = "/"), opts$ntree, opts$maxiter, opts$seed, opts$backend,
                      sep = "|")
  JOB_CONTEXT <<- list(options = opts, jobs = jobs, hashes = hashes, signature = signature)
  WELL_CACHE <<- new.env(parent = emptyenv())
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_table(jobs, file.path(opts$output_dir, "run_manifest.csv"))
  writeLines(c(paste(names(opts), vapply(opts, function(x) paste(x, collapse = ","), character(1)), sep = " = "),
               paste(names(versions), versions, sep = " = ")), file.path(opts$output_dir, "run_settings.txt"))
  writeLines(capture.output(sessionInfo()), file.path(opts$output_dir, "sessionInfo.txt"))
  cat(sprintf("Running %d well/configuration pairs with %d worker(s), %d trees and maxiter=%d.\n",
              nrow(jobs), opts$workers, opts$ntree, opts$maxiter))
  if (opts$workers == 1L) {
    results <- lapply(seq_len(nrow(jobs)), run_one)
  } else {
    cluster <- parallel::makePSOCKcluster(min(opts$workers, nrow(jobs)),
                                         outfile = file.path(opts$output_dir, "workers.log"))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(path) source(path, encoding = "UTF-8"), script)
    parallel::clusterExport(cluster, "JOB_CONTEXT", envir = .GlobalEnv)
    parallel::clusterEvalQ(cluster, {
      data.table::setDTthreads(1L)
      WELL_CACHE <- new.env(parent = emptyenv())
      NULL
    })
    results <- parallel::parLapplyLB(cluster, seq_len(nrow(jobs)), run_one)
  }
  metrics <- do.call(rbind, results)
  metrics <- metrics[order(match(metrics$Configuration, configs$configuration), metrics$Target_Well), ]
  summarize_metrics(metrics, jobs, opts$output_dir)
  cat(sprintf("Finished: %d/%d complete; metrics saved to %s\n",
              sum(metrics$Status == "complete"), nrow(metrics), opts$output_dir))
  if (any(metrics$Status != "complete")) stop("Some jobs did not complete; inspect the saved metrics and resume.")
  invisible(metrics)
}

if (sys.nframe() == 0L) main()
