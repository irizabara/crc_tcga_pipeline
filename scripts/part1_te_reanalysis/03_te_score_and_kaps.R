#!/usr/bin/env Rscript

# ==============================================================================
# TE score generation and KAPS survival partitioning
# ==============================================================================
#
# Expected inputs
# ---------------
# 1. An RDS file containing a numeric TE-expression matrix:
#       rows    = TE subfamilies
#       columns = samples
#    Expression should already be normalized as described in the thesis
#    (RLE-normalized log2 CPM with prior count 5).
#
# 2. An RDS file containing the matched clinical data frame:
#       rows = samples in the same order or identifiable by sample ID
#    Required fields:
#       OS.time (or one of the accepted aliases below)
#       OS      (0/1 event indicator; accepted aliases below)
#       seq_type for the Zhu-score stratified analysis
#
# Usage
# -----
# Rscript 03_te_scores_and_kaps.R \
#   path/to/te_matrix.RDS \
#   path/to/clinical.RDS \
#   path/to/output_directory
#
# The source data are not distributed with this repository.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

if (!requireNamespace("kaps", quietly = TRUE)) {
  stop(
    paste(
      "The 'kaps' package is required.",
      "Install it with:",
      'remotes::install_github("sooheang/kaps", ref = "v1.1.5")',
      sep = "\n"
    )
  )
}

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    paste(
      "Usage:",
      "Rscript 03_te_scores_and_kaps.R <te_matrix.RDS> <clinical.RDS> [output_dir]"
    )
  )
}

te_matrix_path <- args[[1]]
clinical_path <- args[[2]]
output_dir <- if (length(args) >= 3L) args[[3]] else "outputs/part1_kaps"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
n_permutations <- 100L

seven_te_signature <- c(
  "LTR106_Mam",
  "LTR19-int",
  "LTR80A",
  "MER44B",
  "MER57E1",
  "MER65C",
  "Tigger11a"
)

zhu_nine_te_signature <- c(
  "AluSq",
  "HERV1_LTRd",
  "LTR21B",
  "MER57F",
  "MER65C",
  "MER92-int",
  "SVA_C",
  "SVA_F",
  "Tigger12A"
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

first_existing_column <- function(data, candidates, label) {
  hit <- candidates[candidates %in% colnames(data)]

  if (length(hit) == 0L) {
    stop(
      label,
      " column not found. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit[[1]]
}

normalise_te_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

resolve_te_rows <- function(requested_names, available_names) {
  if (is.null(available_names)) {
    stop("The TE-expression matrix must have row names.")
  }

  available_normalised <- normalise_te_name(available_names)
  requested_normalised <- normalise_te_name(requested_names)

  resolved <- vapply(
    seq_along(requested_names),
    function(i) {
      matches <- which(available_normalised == requested_normalised[[i]])

      if (length(matches) == 0L) {
        stop("TE not found in expression matrix: ", requested_names[[i]])
      }

      if (length(matches) > 1L) {
        stop(
          "Ambiguous TE name after normalization: ",
          requested_names[[i]],
          ". Matching rows: ",
          paste(available_names[matches], collapse = ", ")
        )
      }

      available_names[matches]
    },
    character(1)
  )

  data.frame(
    requested_name = requested_names,
    matched_matrix_row = unname(resolved),
    stringsAsFactors = FALSE
  )
}

derive_sample_ids <- function(clinical) {
  direct_candidates <- c(
    "sample_id",
    "sample",
    "Sample",
    "submitter_id",
    "barcode",
    "bcr_sample_barcode"
  )

  direct_hit <- direct_candidates[direct_candidates %in% colnames(clinical)]

  if (length(direct_hit) > 0L) {
    return(as.character(clinical[[direct_hit[[1]]]]))
  }

  if ("bcr_patient_barcode" %in% colnames(clinical)) {
    return(paste0(as.character(clinical$bcr_patient_barcode), "-01A"))
  }

  if (!is.null(rownames(clinical)) &&
      !all(rownames(clinical) == as.character(seq_len(nrow(clinical))))) {
    return(rownames(clinical))
  }

  stop(
    paste(
      "Could not determine clinical sample IDs.",
      "Add a sample_id column or row names matching the TE matrix columns."
    )
  )
}

align_clinical_to_matrix <- function(te_matrix, clinical) {
  matrix_ids <- colnames(te_matrix)

  if (is.null(matrix_ids)) {
    stop("The TE-expression matrix must have column names.")
  }

  clinical_ids <- derive_sample_ids(clinical)
  matched_rows <- match(matrix_ids, clinical_ids)

  if (anyNA(matched_rows)) {
    missing_ids <- matrix_ids[is.na(matched_rows)]

    stop(
      "Clinical data are missing ",
      length(missing_ids),
      " matrix samples. First missing IDs: ",
      paste(head(missing_ids, 10L), collapse = ", ")
    )
  }

  aligned <- clinical[matched_rows, , drop = FALSE]
  aligned$sample_id <- matrix_ids

  aligned
}

generate_mean_te_score <- function(te_matrix, signature) {
  mapping <- resolve_te_rows(signature, rownames(te_matrix))
  signature_matrix <- te_matrix[mapping$matched_matrix_row, , drop = FALSE]

  if (!all(vapply(signature_matrix, is.numeric, logical(1)))) {
    storage.mode(signature_matrix) <- "numeric"
  }

  score <- colMeans(signature_matrix, na.rm = TRUE)

  if (any(!is.finite(score))) {
    stop("Non-finite TE scores were produced.")
  }

  list(
    score = score,
    mapping = mapping
  )
}

prepare_survival_data <- function(clinical, score, score_name) {
  os_time_col <- first_existing_column(
    clinical,
    c("OS.time", "OS_time", "os_time", "overall_survival_time", "time"),
    "Overall-survival time"
  )

  os_event_col <- first_existing_column(
    clinical,
    c("OS", "OS_event", "os_event", "overall_survival_event", "status"),
    "Overall-survival event"
  )

  data <- data.frame(
    sample_id = clinical$sample_id,
    os_time = as.numeric(clinical[[os_time_col]]),
    os_event = as.numeric(as.character(clinical[[os_event_col]])),
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )

  if ("seq_type" %in% colnames(clinical)) {
    data$seq_type <- as.character(clinical$seq_type)
  }

  data <- data[complete.cases(data[, c("os_time", "os_event", "score")]), ]

  if (!all(data$os_event %in% c(0, 1))) {
    stop(score_name, ": OS event values must be encoded as 0/1.")
  }

  if (nrow(data) < 20L || sum(data$os_event) < 5L) {
    stop(score_name, ": insufficient complete survival data.")
  }

  data
}

extract_candidate_fit <- function(search_fit, candidate_k, tested_k) {
  if (length(tested_k) == 1L) {
    return(search_fit)
  }

  candidate_index <- match(candidate_k, tested_k)

  if (is.na(candidate_index)) {
    stop("Requested K was not included in the KAPS search.")
  }

  search_fit@results[[candidate_index]]
}

group_median_survival <- function(survival_fit) {
  fit_table <- summary(survival_fit)$table

  if (is.null(dim(fit_table))) {
    fit_table <- matrix(
      fit_table,
      nrow = 1L,
      dimnames = list(names(survival_fit$strata), names(fit_table))
    )
  }

  result <- data.frame(
    stratum = rownames(fit_table),
    median_survival = as.numeric(fit_table[, "median"]),
    stringsAsFactors = FALSE
  )

  result$group <- sub("^group=", "", result$stratum)
  result[, c("group", "median_survival")]
}

summarise_candidate_k <- function(candidate_fit, analysis_data, analysis, cohort, k) {
  group_id <- as.integer(candidate_fit@groupID)

  if (length(group_id) != nrow(analysis_data)) {
    predicted <- predict(candidate_fit, newdata = analysis_data, type = "predict")
    group_id <- as.integer(predicted$Group)
  }

  ordered_groups <- sort(unique(group_id))
  group_labels <- paste0("G", seq_along(ordered_groups))
  names(group_labels) <- as.character(ordered_groups)

  labelled_data <- analysis_data
  labelled_data$group <- factor(
    unname(group_labels[as.character(group_id)]),
    levels = paste0("G", seq_len(k))
  )

  logrank_fit <- survival::survdiff(
    survival::Surv(os_time, os_event) ~ group,
    data = labelled_data
  )

  logrank_df <- max(length(logrank_fit$n) - 1L, 1L)
  logrank_p <- stats::pchisq(
    logrank_fit$chisq,
    df = logrank_df,
    lower.tail = FALSE
  )

  km_fit <- survival::survfit(
    survival::Surv(os_time, os_event) ~ group,
    data = labelled_data
  )

  median_table <- group_median_survival(km_fit)

  group_summary <- labelled_data |>
    group_by(group) |>
    summarise(
      n = n(),
      events = sum(os_event),
      score_min = min(score),
      score_median = median(score),
      score_max = max(score),
      .groups = "drop"
    ) |>
    left_join(median_table, by = "group") |>
    mutate(
      analysis = analysis,
      cohort = cohort,
      K = k,
      .before = 1
    )

  cut_points <- as.numeric(candidate_fit@split.pt)

  overall_summary <- data.frame(
    analysis = analysis,
    cohort = cohort,
    K = k,
    n = nrow(labelled_data),
    events = sum(labelled_data$os_event),
    logrank_chisq = unname(logrank_fit$chisq),
    logrank_df = logrank_df,
    logrank_p = logrank_p,
    cut_points = paste(signif(cut_points, 8L), collapse = ";"),
    stringsAsFactors = FALSE
  )

  assignments <- labelled_data |>
    mutate(
      analysis = analysis,
      cohort = cohort,
      K = k,
      .before = 1
    )

  list(
    overall = overall_summary,
    groups = group_summary,
    assignments = assignments,
    km_fit = km_fit
  )
}

plot_km <- function(km_fit, path, title) {
  grDevices::pdf(path, width = 7, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  n_groups <- length(km_fit$strata)

  plot(
    km_fit,
    lty = seq_len(n_groups),
    lwd = 2,
    xlab = "Overall survival time (days)",
    ylab = "Survival probability",
    main = title,
    mark.time = TRUE
  )

  legend(
    "bottomleft",
    legend = names(km_fit$strata),
    lty = seq_len(n_groups),
    lwd = 2,
    bty = "n"
  )
}

run_kaps_search <- function(
    analysis_data,
    analysis,
    cohort,
    k_values,
    n_permutations,
    output_dir,
    seed
) {
  analysis_data <- analysis_data |>
    arrange(score)

  set.seed(seed)

  search_fit <- kaps::kaps(
    survival::Surv(os_time, os_event) ~ score,
    data = analysis_data,
    K = k_values,
    type = "perm",
    N.perm = n_permutations
  )

  saveRDS(
    search_fit,
    file.path(
      output_dir,
      paste0(analysis, "_", cohort, "_kaps_fit_K", min(k_values), "-", max(k_values), ".RDS")
    )
  )

  overall_results <- list()
  group_results <- list()

  for (k in k_values) {
    candidate_fit <- extract_candidate_fit(search_fit, k, k_values)

    candidate_summary <- summarise_candidate_k(
      candidate_fit = candidate_fit,
      analysis_data = analysis_data,
      analysis = analysis,
      cohort = cohort,
      k = k
    )

    overall_results[[as.character(k)]] <- candidate_summary$overall
    group_results[[as.character(k)]] <- candidate_summary$groups

    write.csv(
      candidate_summary$assignments,
      file.path(
        output_dir,
        paste0(analysis, "_", cohort, "_K", k, "_assignments.csv")
      ),
      row.names = FALSE
    )

    plot_km(
      candidate_summary$km_fit,
      file.path(
        output_dir,
        paste0(analysis, "_", cohort, "_K", k, "_KM.pdf")
      ),
      paste0(analysis, " — ", cohort, " — K=", k)
    )
  }

  list(
    fit = search_fit,
    overall = bind_rows(overall_results),
    groups = bind_rows(group_results)
  )
}

# ------------------------------------------------------------------------------
# Load and align inputs
# ------------------------------------------------------------------------------

te_matrix <- readRDS(te_matrix_path)
clinical <- readRDS(clinical_path)

if (!is.matrix(te_matrix) && !is.data.frame(te_matrix)) {
  stop("The TE-expression RDS must contain a matrix or data frame.")
}

te_matrix <- as.matrix(te_matrix)
storage.mode(te_matrix) <- "numeric"

clinical <- align_clinical_to_matrix(te_matrix, clinical)

cat("Aligned samples:", nrow(clinical), "\n")
cat("Available TE subfamilies:", nrow(te_matrix), "\n")

# ------------------------------------------------------------------------------
# Analysis 1: Confounder-aware seven-TE score, full cohort, K = 2 to 5
# ------------------------------------------------------------------------------

seven_score_result <- generate_mean_te_score(
  te_matrix = te_matrix,
  signature = seven_te_signature
)

write.csv(
  seven_score_result$mapping,
  file.path(output_dir, "seven_te_signature_mapping.csv"),
  row.names = FALSE
)

seven_data <- prepare_survival_data(
  clinical = clinical,
  score = seven_score_result$score,
  score_name = "Seven-TE score"
)

write.csv(
  seven_data,
  file.path(output_dir, "seven_te_scores.csv"),
  row.names = FALSE
)

seven_kaps <- run_kaps_search(
  analysis_data = seven_data,
  analysis = "seven_TE_score",
  cohort = "ALL",
  k_values = 2:5,
  n_permutations = n_permutations,
  output_dir = output_dir,
  seed = 2026
)

# ------------------------------------------------------------------------------
# Analysis 2: Zhu nine-TE score and sequencing-type-stratified KAPS
# ------------------------------------------------------------------------------

zhu_score_result <- generate_mean_te_score(
  te_matrix = te_matrix,
  signature = zhu_nine_te_signature
)

write.csv(
  zhu_score_result$mapping,
  file.path(output_dir, "zhu_nine_te_signature_mapping.csv"),
  row.names = FALSE
)

zhu_data <- prepare_survival_data(
  clinical = clinical,
  score = zhu_score_result$score,
  score_name = "Zhu nine-TE score"
)

if (!"seq_type" %in% colnames(zhu_data)) {
  stop(
    paste(
      "The clinical data must contain seq_type to reproduce",
      "the sequencing-stratified Zhu-score KAPS analysis."
    )
  )
}

write.csv(
  zhu_data,
  file.path(output_dir, "zhu_nine_te_scores.csv"),
  row.names = FALSE
)

cohort_definitions <- list(
  ALL = rep(TRUE, nrow(zhu_data)),
  PE = toupper(zhu_data$seq_type) %in% c("PE", "PAIRED", "PAIRED-END", "PAIRED_END"),
  SE = toupper(zhu_data$seq_type) %in% c("SE", "SINGLE", "SINGLE-END", "SINGLE_END")
)

zhu_overall_results <- list()
zhu_group_results <- list()

for (cohort_name in names(cohort_definitions)) {
  cohort_data <- zhu_data[cohort_definitions[[cohort_name]], , drop = FALSE]

  if (nrow(cohort_data) == 0L) {
    warning("Skipping empty cohort: ", cohort_name)
    next
  }

  cohort_result <- run_kaps_search(
    analysis_data = cohort_data,
    analysis = "zhu_nine_TE_score",
    cohort = cohort_name,
    k_values = 2:4,
    n_permutations = n_permutations,
    output_dir = output_dir,
    seed = 2026 + match(cohort_name, names(cohort_definitions))
  )

  zhu_overall_results[[cohort_name]] <- cohort_result$overall
  zhu_group_results[[cohort_name]] <- cohort_result$groups
}

# ------------------------------------------------------------------------------
# Save aggregate summaries
# ------------------------------------------------------------------------------

all_overall_results <- bind_rows(
  seven_kaps$overall,
  bind_rows(zhu_overall_results)
)

all_group_results <- bind_rows(
  seven_kaps$groups,
  bind_rows(zhu_group_results)
)

write.csv(
  all_overall_results,
  file.path(output_dir, "kaps_overall_summary.csv"),
  row.names = FALSE
)

write.csv(
  all_group_results,
  file.path(output_dir, "kaps_group_summary.csv"),
  row.names = FALSE
)

cat("\nKAPS analyses complete.\n")
cat("Outputs written to:", normalizePath(output_dir), "\n")
cat("\nOverall summary:\n")
print(all_overall_results)
