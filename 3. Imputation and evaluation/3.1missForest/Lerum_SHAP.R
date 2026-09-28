options_from_args <- function(args, here) {
  opts <- list(well = "Lerum_1", input_dir = file.path(dirname(dirname(here)),
    "2.data_preprocessing", "2.2 missForest", "Lerum"),
    output_dir = file.path(here, "Lerum_SHAP"), max_rows = 3000,
    nsim = 30L, only_removed = FALSE)
  if ("--help" %in% args) {
    cat("Rscript np_shap_lerum.R [--well Lerum_1] [--max-rows 3000|Inf]\n",
        "  [--nsim 30] [--only-removed] [--input-dir PATH] [--output-dir PATH]\n")
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    if (args[i] == "--only-removed") {
      opts$only_removed <- TRUE
      i <- i + 1L
      next
    }
    key <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
    if (!startsWith(args[i], "--") || !key %in% setdiff(names(opts), "only_removed") ||
        i == length(args)) stop("Invalid option: ", args[i])
    opts[[key]] <- args[i + 1L]
    i <- i + 2L
  }
  opts$max_rows <- suppressWarnings(as.numeric(opts$max_rows))
  opts$nsim <- suppressWarnings(as.numeric(opts$nsim))
  if (is.na(opts$max_rows) || opts$max_rows < 2 || opts$max_rows != floor(opts$max_rows) ||
      !is.finite(opts$nsim) || opts$nsim < 2 || opts$nsim != floor(opts$nsim)) {
    stop("max-rows must be an integer >= 2 or Inf; nsim must be an integer >= 2.")
  }
  if (!grepl("^[A-Za-z0-9_]+$", opts$well)) stop("Invalid well identifier.")
  opts
}

read_csv <- function(path) read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
predict_rf <- function(object, newdata) predict(object, data = newdata, num.threads = 1L)$predictions

main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), winslash = "/")
  opts <- options_from_args(commandArgs(trailingOnly = TRUE), dirname(script))
  if (is.null(opts)) return(invisible(NULL))
  packages <- c("ranger", "fastshap", "shapviz", "ggplot2")
  for (p in packages) if (!requireNamespace(p, quietly = TRUE)) stop("Install required R package: ", p)
  suppressPackageStartupMessages(library(ggplot2))
  features <- paste0("NP_", 1:12, "m")
  input_file <- file.path(opts$input_dir, "data", paste0(opts$well, ".csv"))
  mask_file <- file.path(opts$input_dir, "mask_records", paste0(opts$well, ".csv"))
  data <- read_csv(input_file)
  mask <- read_csv(mask_file)
  if (length(setdiff(c("date", "target_gwl_original", features), names(data))) ||
      length(setdiff(c("date", "was_removed"), names(mask)))) stop("Missing input columns.")
  data$date <- as.Date(data$date)
  mask$date <- as.Date(mask$date)
  if (anyNA(data$date) || anyDuplicated(data$date) || !identical(data$date, mask$date)) {
    stop("Data and removal records must contain the same unique dates in the same order.")
  }
  flags <- tolower(trimws(as.character(mask$was_removed)))
  if (anyNA(flags) || any(!flags %in% c("true", "false", "1", "0"))) stop("Invalid removal flags.")
  data$was_removed <- flags %in% c("true", "1")
  data <- data[order(data$date), ]
  if (!all(vapply(data[c("target_gwl_original", features)], is.numeric, logical(1)))) {
    stop("Target and NP predictors must be numeric.")
  }
  if (any(data$was_removed & !is.na(data$target_gwl_original))) stop("Held-out target values are not masked.")
  if (any(is.infinite(as.matrix(data[c("target_gwl_original", features)])))) stop("Infinite model inputs.")
  complete_np <- complete.cases(data[features])
  training <- !data$was_removed & is.finite(data$target_gwl_original) & complete_np
  if (sum(training) < 100 || length(unique(data$target_gwl_original[training])) < 2) {
    stop("Need at least 100 complete observed training rows and a varying target.")
  }
  indices <- which(complete_np & (!opts$only_removed | data$was_removed))
  eligible <- length(indices)
  if (!eligible) stop("No dates with complete NP predictors to explain.")
  if (eligible > opts$max_rows) {
    indices <- indices[unique(round(seq(1, eligible, length.out = opts$max_rows)))]
  }
  background <- data[training, features, drop = FALSE]
  explain <- data[indices, features, drop = FALSE]
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(42)

  model <- ranger::ranger(x = background, y = data$target_gwl_original[training],
    num.trees = 300L, mtry = 3L, min.node.size = 5L,
    importance = "permutation", seed = 42L, num.threads = 1L)
  predictions <- predict_rf(model, explain)
  baseline <- mean(predict_rf(model, background))
  cat(sprintf("%s: %d training rows; explaining %d/%d eligible dates, nsim=%d.\n",
              opts$well, sum(training), length(indices), eligible, opts$nsim))
  set.seed(42)
  shap <- fastshap::explain(model, X = background, newdata = explain,
    pred_wrapper = predict_rf, nsim = opts$nsim, adjust = TRUE, baseline = baseline)
  shap <- as.matrix(shap)[, features, drop = FALSE]
  residual <- predictions - baseline - rowSums(shap)
  tolerance <- 1e-7 * max(1, max(abs(predictions - baseline)))
  if (any(!is.finite(shap)) || max(abs(residual)) > tolerance) {
    stop("SHAP values are non-finite or fail prediction = baseline + sum(SHAP).")
  }
  output <- file.path(opts$output_dir, opts$well)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  wide <- data.frame(date = data$date[indices], was_removed = data$was_removed[indices],
    prediction = predictions, baseline = baseline, shap, Total_NP_SHAP = rowSums(shap),
    additivity_error = residual, check.names = FALSE)
  ranking <- data.frame(NP_Period = features, Mean_Absolute_SHAP = colMeans(abs(shap)))
  ranking <- ranking[order(-ranking$Mean_Absolute_SHAP), ]
  ranking$Rank <- seq_len(nrow(ranking))
  write_csv(wide, file.path(output, "NP_SHAP_values.csv"))
  write_csv(data.frame(date = data$date[indices], explain), file.path(output, "NP_explained_features.csv"))
  write_csv(ranking, file.path(output, "NP_SHAP_importance_ranking.csv"))
  write_csv(data.frame(date = data$date, used_for_training = training,
    eligible_for_explanation = complete_np & (!opts$only_removed | data$was_removed),
    selected_for_explanation = seq_len(nrow(data)) %in% indices), file.path(output, "row_selection.csv"))

  scope <- if (opts$only_removed) "artificial gaps" else "the full eligible record"
  subtitle <- sprintf("Separate NP random forest | %d dates sampled from %s", length(indices), scope)
  ranking$NP_Period <- factor(ranking$NP_Period, levels = rev(ranking$NP_Period))
  importance_plot <- ggplot(ranking, aes(Mean_Absolute_SHAP, NP_Period)) +
    geom_col(fill = "#317A9B", width = 0.7) + theme_minimal(base_size = 11) +
    labs(title = paste(opts$well, "- NP-period importance"), subtitle = subtitle,
         x = "Mean absolute SHAP value", y = "NP rolling period")
  long <- data.frame(date = rep(data$date[indices], times = 12),
    NP_Period = factor(rep(features, each = length(indices)), levels = features),
    SHAP_Value = as.vector(shap))

  layer <- if (opts$only_removed) geom_point(size = 0.7, alpha = 0.65) else geom_line(linewidth = 0.3, alpha = 0.7)
  time_plot <- ggplot(long, aes(date, SHAP_Value, color = NP_Period)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") + layer +
    theme_minimal(base_size = 11) + labs(title = paste(opts$well, "- NP SHAP over time"),
      subtitle = subtitle, x = "Date", y = "SHAP contribution to predicted water level", color = "NP period",
      caption = "Black: sum of all 12 NP contributions relative to the background prediction.")
  if (opts$only_removed) {
    time_plot <- time_plot + geom_point(data = wide, aes(date, Total_NP_SHAP),
                                        inherit.aes = FALSE, size = 0.7)
  } else {
    time_plot <- time_plot + geom_line(data = wide, aes(date, Total_NP_SHAP),
                                       inherit.aes = FALSE, linewidth = 0.55)
  }
  shap_object <- shapviz::shapviz(shap, X = explain, baseline = baseline)
  set.seed(42)
  beeswarm_plot <- shapviz::sv_importance(shap_object, kind = "beeswarm", max_display = 12) +
    labs(title = paste(opts$well, "- NP SHAP summary"), subtitle = subtitle) +
    theme(plot.title = element_text(size = 13), plot.subtitle = element_text(size = 10))
  ggsave(file.path(output, "NP_SHAP_importance.png"), importance_plot, width = 10, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(output, "NP_SHAP_over_time.png"), time_plot, width = 14, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(output, "NP_SHAP_beeswarm.png"), beeswarm_plot, width = 10, height = 6, dpi = 300, bg = "white")
  provenance <- c(unlist(opts), num_trees = 300, mtry = 3, min_node_size = 5, seed = 42,
    num_threads = 1, training_rows = sum(training), eligible_rows = eligible,
    explained_rows = length(indices), baseline = baseline, OOB_MSE = model$prediction.error,
    max_additivity_error = max(abs(residual)), additivity_tolerance = tolerance,
    input_md5 = unname(tools::md5sum(input_file)), mask_md5 = unname(tools::md5sum(mask_file)),
    script_md5 = unname(tools::md5sum(script)),
    setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), packages))
  writeLines(paste(names(provenance), provenance, sep = " = "), file.path(output, "run_settings.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  cat("Saved SHAP tables and three plots to", normalizePath(output, winslash = "/"), "\n")
}

if (sys.nframe() == 0L) main()
