#!/usr/bin/env Rscript

# ==============================================================================
# Part 3 skeleton: MSS paired-end TE-PC1 analysis
# ==============================================================================
#
# Purpose
# -------
# This script documents the main analytical logic used in Part 3 of the thesis:
#
#   1. restrict the cohort to microsatellite-stable (MSS), paired-end samples;
#   2. perform PCA on normalized TE expression;
#   3. interpret PC1 through TE-family loadings;
#   4. test PC1 associations with survival and immune/stromal indices;
#   5. generate median-split and extreme-quartile labels for CLAM;
#   6. create a host-gene ranking for GSEA;
#   7. optionally join HoVer-Net cell-composition summaries.
#
# This is a curated workflow template. It is designed to make the thesis logic
# inspectable rather than reproduce institution-specific paths or every
# exploratory analysis exactly.
#
# Expected inputs
# ---------------
# Required:
#   - clinical RDS: one row per patient/sample
#   - TE-expression RDS: rows = TE subfamilies, columns = samples
#
# Optional:
#   - immune/stromal score RDS: rows = samples, columns = scores
#   - host-gene expression RDS: rows = genes, columns = samples
#   - slide metadata CSV: slide_id and patient_id columns
#   - HoVer-Net composition CSV produced by the Part 2 workflow
#
# Expression matrices should already contain the normalized values used in the
# thesis. Raw read processing and REdiscoverTE quantification are not repeated
# here.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

# ------------------------------------------------------------------------------
# Configuration: replace placeholders with local paths
# ------------------------------------------------------------------------------

clinical_path <- "<path_to_clinical_data.RDS>"
te_matrix_path <- "<path_to_TE_expression_matrix.RDS>"

# Optional inputs. Leave as NA_character_ to skip the relevant section.
immune_scores_path <- NA_character_
gene_expression_path <- NA_character_
slide_metadata_path <- NA_character_
hovernet_composition_path <- NA_character_

output_dir <- "outputs/part3_mss_pe_pc1"

# Column-name candidates are resolved automatically where possible.
preferred_pc1_name <- "PC1"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------------------

first_existing_column <- function(data, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% colnames(data)]

  if (length(hit) == 0L) {
    if (required) {
      stop(label, " column not found. Tried: ", paste(candidates, collapse = ", "))
    }
    return(NULL)
  }

  hit[[1]]
}

derive_sample_ids <- function(data) {
  candidates <- c(
    "sample_id",
    "sample",
    "Sample",
    "submitter_id",
    "barcode",
    "bcr_sample_barcode"
  )

  hit <- candidates[candidates %in% colnames(data)]

  if (length(hit) > 0L) {
    return(as.character(data[[hit[[1]]]]))
  }

  if ("bcr_patient_barcode" %in% colnames(data)) {
    return(paste0(as.character(data$bcr_patient_barcode), "-01A"))
  }

  if (!is.null(rownames(data)) &&
      !all(rownames(data) == as.character(seq_len(nrow(data))))) {
    return(rownames(data))
  }

  stop(
    paste(
      "Could not determine sample IDs.",
      "Provide a sample_id column or informative row names."
    )
  )
}

derive_patient_ids <- function(data, sample_ids = NULL) {
  candidates <- c(
    "patient_id",
    "bcr_patient_barcode",
    "case_submitter_id",
    "case_id"
  )

  hit <- candidates[candidates %in% colnames(data)]

  if (length(hit) > 0L) {
    return(as.character(data[[hit[[1]]]]))
  }

  if (is.null(sample_ids)) {
    sample_ids <- derive_sample_ids(data)
  }

  # TCGA patient barcode: TCGA-XX-YYYY
  sub("^(([^-]+-){2}[^-]+).*", "\\1", sample_ids)
}

clean_stage <- function(x) {
  x <- trimws(as.character(x))

  stage <- case_when(
    grepl("^Stage IV|^IV", x, ignore.case = TRUE) ~ "IV",
    grepl("^Stage III|^III", x, ignore.case = TRUE) ~ "III",
    grepl("^Stage II|^II", x, ignore.case = TRUE) ~ "II",
    grepl("^Stage I|^I", x, ignore.case = TRUE) ~ "I",
    TRUE ~ NA_character_
  )

  factor(stage, levels = c("I", "II", "III", "IV"))
}

coerce_event <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  output <- rep(NA_real_, length(x_chr))

  output[x_chr %in% c("1", "true", "yes", "event", "dead", "deceased")] <- 1
  output[x_chr %in% c("0", "false", "no", "censored", "alive")] <- 0

  numeric_values <- suppressWarnings(as.numeric(x_chr))
  numeric_ok <- is.na(output) & numeric_values %in% c(0, 1)
  output[numeric_ok] <- numeric_values[numeric_ok]

  output
}

path_is_available <- function(path) {
  !is.na(path) && nzchar(path) && file.exists(path)
}

# ------------------------------------------------------------------------------
# Cohort filtering and matrix alignment
# ------------------------------------------------------------------------------

is_mss_value <- function(x) {
  normalized <- toupper(trimws(as.character(x)))

  normalized %in% c(
    "MSS",
    "MICROSATELLITE STABLE",
    "MICROSATELLITE-STABLE",
    "STABLE"
  )
}

is_paired_end_value <- function(x) {
  normalized <- toupper(trimws(as.character(x)))

  normalized %in% c(
    "PE",
    "PAIRED",
    "PAIRED-END",
    "PAIRED_END",
    "PAIRED END"
  )
}

align_clinical_to_matrix <- function(clinical, matrix) {
  if (is.null(colnames(matrix))) {
    stop("Expression matrix must have sample IDs as column names.")
  }

  clinical_ids <- derive_sample_ids(clinical)
  matrix_ids <- colnames(matrix)
  matched_rows <- match(matrix_ids, clinical_ids)

  if (anyNA(matched_rows)) {
    missing_ids <- matrix_ids[is.na(matched_rows)]

    stop(
      "Clinical data are missing ",
      length(missing_ids),
      " expression-matrix samples. First missing IDs: ",
      paste(head(missing_ids, 10L), collapse = ", ")
    )
  }

  aligned <- clinical[matched_rows, , drop = FALSE]
  aligned$sample_id <- matrix_ids
  aligned$patient_id <- derive_patient_ids(aligned, matrix_ids)
  aligned
}

subset_matrix_to_samples <- function(matrix, sample_ids) {
  missing_ids <- setdiff(sample_ids, colnames(matrix))

  if (length(missing_ids) > 0L) {
    stop(
      "Matrix is missing ",
      length(missing_ids),
      " requested samples. First missing IDs: ",
      paste(head(missing_ids, 10L), collapse = ", ")
    )
  }

  matrix[, sample_ids, drop = FALSE]
}

# ------------------------------------------------------------------------------
# PCA and TE-family interpretation
# ------------------------------------------------------------------------------

impute_feature_medians <- function(sample_by_feature) {
  for (column_index in seq_len(ncol(sample_by_feature))) {
    values <- sample_by_feature[, column_index]

    if (anyNA(values)) {
      feature_median <- median(values, na.rm = TRUE)

      if (!is.finite(feature_median)) {
        stop(
          "A TE feature contains only missing or non-finite values: ",
          colnames(sample_by_feature)[column_index]
        )
      }

      values[is.na(values)] <- feature_median
      sample_by_feature[, column_index] <- values
    }
  }

  sample_by_feature
}

classify_te_family <- function(te_name) {
  name_upper <- toupper(te_name)

  case_when(
    grepl("ALU|SINE", name_upper) ~ "SINE/Alu",
    grepl("ERV|HERV|LTR", name_upper) ~ "LTR/ERV",
    grepl("LINE|(^|[^A-Z])L1([^A-Z]|$)|(^|[^A-Z])L2([^A-Z]|$)", name_upper) ~
      "LINE/L1/L2",
    grepl("MER", name_upper) ~ "MER",
    grepl("TIGGER|DNA", name_upper) ~ "DNA transposon",
    TRUE ~ "Other"
  )
}

orient_pc1_for_interpretation <- function(scores, loadings) {
  # PCA signs are arbitrary. For consistent reporting, orient PC1 so that the
  # positive direction corresponds more strongly to the SINE/Alu side than to
  # the combined LTR/ERV and LINE/L1/L2 side, matching the thesis interpretation.

  families <- classify_te_family(rownames(loadings))
  pc1_loadings <- loadings[, 1]

  sine_mean <- mean(pc1_loadings[families == "SINE/Alu"], na.rm = TRUE)
  retro_mean <- mean(
    pc1_loadings[families %in% c("LTR/ERV", "LINE/L1/L2")],
    na.rm = TRUE
  )

  should_flip <- is.finite(sine_mean) &&
    is.finite(retro_mean) &&
    sine_mean < retro_mean

  if (should_flip) {
    scores[, 1] <- -scores[, 1]
    loadings[, 1] <- -loadings[, 1]
  }

  list(
    scores = scores,
    loadings = loadings,
    flipped = should_flip
  )
}

run_te_pca <- function(te_matrix) {
  finite_variance <- apply(
    te_matrix,
    1,
    function(values) {
      variance <- var(as.numeric(values), na.rm = TRUE)
      is.finite(variance) && variance > 0
    }
  )

  filtered <- te_matrix[finite_variance, , drop = FALSE]
  sample_by_te <- t(filtered)
  sample_by_te <- impute_feature_medians(sample_by_te)

  pca <- prcomp(
    sample_by_te,
    center = TRUE,
    scale. = TRUE
  )

  oriented <- orient_pc1_for_interpretation(
    scores = pca$x,
    loadings = pca$rotation
  )

  pca$x <- oriented$scores
  pca$rotation <- oriented$loadings

  variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  list(
    pca = pca,
    variance_explained = variance_explained,
    pc1_flipped = oriented$flipped
  )
}

summarize_pc1_loadings <- function(loadings) {
  loading_table <- data.frame(
    TE = rownames(loadings),
    PC1_loading = as.numeric(loadings[, 1]),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      family = classify_te_family(TE),
      loading_direction = case_when(
        PC1_loading > 0 ~ "positive_PC1",
        PC1_loading < 0 ~ "negative_PC1",
        TRUE ~ "zero"
      ),
      absolute_loading = abs(PC1_loading)
    ) |>
    arrange(desc(absolute_loading))

  family_summary <- loading_table |>
    group_by(family, loading_direction) |>
    summarise(
      n_TEs = n(),
      mean_loading = mean(PC1_loading),
      median_loading = median(PC1_loading),
      mean_absolute_loading = mean(absolute_loading),
      .groups = "drop"
    )

  list(
    loading_table = loading_table,
    family_summary = family_summary
  )
}

# ------------------------------------------------------------------------------
# Survival analysis
# ------------------------------------------------------------------------------

run_pc1_survival <- function(analysis_table) {
  os_time_col <- first_existing_column(
    analysis_table,
    c("OS.time", "OS_time", "os_time"),
    "Overall-survival time"
  )

  os_event_col <- first_existing_column(
    analysis_table,
    c("OS", "OS_event", "os_event"),
    "Overall-survival event"
  )

  age_col <- first_existing_column(
    analysis_table,
    c("age", "age_at_initial_pathologic_diagnosis"),
    "Age"
  )

  stage_col <- first_existing_column(
    analysis_table,
    c("stage", "ajcc_pathologic_tumor_stage"),
    "Stage"
  )

  survival_data <- data.frame(
    patient_id = analysis_table$patient_id,
    sample_id = analysis_table$sample_id,
    os_time = suppressWarnings(as.numeric(analysis_table[[os_time_col]])),
    os_event = coerce_event(analysis_table[[os_event_col]]),
    age = suppressWarnings(as.numeric(analysis_table[[age_col]])),
    stage = clean_stage(analysis_table[[stage_col]]),
    PC1 = analysis_table$PC1,
    stringsAsFactors = FALSE
  )

  survival_data <- survival_data[
    complete.cases(survival_data[, c(
      "os_time",
      "os_event",
      "age",
      "stage",
      "PC1"
    )]),
    ,
    drop = FALSE
  ]

  continuous_fit <- coxph(
    Surv(os_time, os_event) ~ PC1 + age + stage,
    data = survival_data
  )

  continuous_summary <- summary(continuous_fit)
  pc1_coefficient <- continuous_summary$coefficients["PC1", , drop = FALSE]
  pc1_confidence <- continuous_summary$conf.int["PC1", , drop = FALSE]

  continuous_result <- data.frame(
    model = "Surv(OS) ~ PC1 + age + stage",
    n = nrow(survival_data),
    events = sum(survival_data$os_event),
    coefficient = unname(pc1_coefficient[1, "coef"]),
    HR = unname(pc1_confidence[1, "exp(coef)"]),
    lower_CI = unname(pc1_confidence[1, "lower .95"]),
    upper_CI = unname(pc1_confidence[1, "upper .95"]),
    p_value = unname(pc1_coefficient[1, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )

  survival_data$pc1_median_group <- factor(
    ifelse(
      survival_data$PC1 >= median(survival_data$PC1),
      "high_PC1",
      "low_PC1"
    ),
    levels = c("low_PC1", "high_PC1")
  )

  logrank_fit <- survdiff(
    Surv(os_time, os_event) ~ pc1_median_group,
    data = survival_data
  )

  logrank_result <- data.frame(
    comparison = "PC1 median split",
    n = nrow(survival_data),
    events = sum(survival_data$os_event),
    chisq = unname(logrank_fit$chisq),
    df = length(logrank_fit$n) - 1L,
    p_value = pchisq(
      logrank_fit$chisq,
      df = length(logrank_fit$n) - 1L,
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )

  list(
    continuous = continuous_result,
    logrank = logrank_result,
    survival_data = survival_data,
    continuous_fit = continuous_fit
  )
}

# ------------------------------------------------------------------------------
# Immune and stromal correlations
# ------------------------------------------------------------------------------

align_sample_table <- function(table, target_sample_ids) {
  table <- as.data.frame(table, check.names = FALSE)
  table_ids <- derive_sample_ids(table)
  matched_rows <- match(target_sample_ids, table_ids)

  if (anyNA(matched_rows)) {
    missing_ids <- target_sample_ids[is.na(matched_rows)]

    stop(
      "Score table is missing ",
      length(missing_ids),
      " cohort samples. First missing IDs: ",
      paste(head(missing_ids, 10L), collapse = ", ")
    )
  }

  aligned <- table[matched_rows, , drop = FALSE]
  aligned$sample_id <- target_sample_ids
  aligned
}

run_pc1_score_correlations <- function(pc1_table, score_table) {
  id_columns <- intersect(
    colnames(score_table),
    c(
      "sample_id",
      "sample",
      "Sample",
      "submitter_id",
      "barcode",
      "bcr_sample_barcode",
      "bcr_patient_barcode",
      "patient_id"
    )
  )

  numeric_columns <- setdiff(
    colnames(score_table)[vapply(score_table, is.numeric, logical(1))],
    id_columns
  )

  results <- lapply(
    numeric_columns,
    function(score_name) {
      values <- score_table[[score_name]]
      complete_index <- complete.cases(pc1_table$PC1, values)

      if (sum(complete_index) < 10L ||
          length(unique(values[complete_index])) < 2L) {
        return(
          data.frame(
            score = score_name,
            n = sum(complete_index),
            rho = NA_real_,
            p_value = NA_real_
          )
        )
      }

      test <- suppressWarnings(
        cor.test(
          pc1_table$PC1[complete_index],
          values[complete_index],
          method = "spearman",
          exact = FALSE
        )
      )

      data.frame(
        score = score_name,
        n = sum(complete_index),
        rho = unname(test$estimate),
        p_value = test$p.value
      )
    }
  )

  bind_rows(results) |>
    mutate(FDR = p.adjust(p_value, method = "BH")) |>
    arrange(p_value)
}

# ------------------------------------------------------------------------------
# CLAM label generation
# ------------------------------------------------------------------------------

make_pc1_groups <- function(pc1_table) {
  median_value <- median(pc1_table$PC1, na.rm = TRUE)
  lower_quartile <- quantile(pc1_table$PC1, 0.25, na.rm = TRUE)
  upper_quartile <- quantile(pc1_table$PC1, 0.75, na.rm = TRUE)

  pc1_table |>
    mutate(
      pc1_median_group = ifelse(
        PC1 >= median_value,
        "high_PC1",
        "low_PC1"
      ),
      pc1_median_label = ifelse(pc1_median_group == "high_PC1", 1L, 0L),
      pc1_extreme_group = case_when(
        PC1 <= lower_quartile ~ "low_PC1",
        PC1 >= upper_quartile ~ "high_PC1",
        TRUE ~ NA_character_
      ),
      pc1_extreme_label = case_when(
        pc1_extreme_group == "low_PC1" ~ 0L,
        pc1_extreme_group == "high_PC1" ~ 1L,
        TRUE ~ NA_integer_
      )
    )
}

prepare_slide_labels <- function(pc1_groups, slide_metadata) {
  slide_id_col <- first_existing_column(
    slide_metadata,
    c("slide_id", "slide", "filename", "slide_filename"),
    "Slide ID"
  )

  patient_id_col <- first_existing_column(
    slide_metadata,
    c("patient_id", "bcr_patient_barcode", "case_submitter_id"),
    "Patient ID"
  )

  slides <- slide_metadata |>
    transmute(
      slide_id = as.character(.data[[slide_id_col]]),
      patient_id = as.character(.data[[patient_id_col]])
    ) |>
    distinct()

  joined <- slides |>
    inner_join(
      pc1_groups |>
        select(
          patient_id,
          PC1,
          pc1_median_group,
          pc1_median_label,
          pc1_extreme_group,
          pc1_extreme_label
        ),
      by = "patient_id"
    )

  median_labels <- joined |>
    transmute(
      slide_id,
      patient_id,
      PC1,
      label = pc1_median_label,
      label_name = pc1_median_group
    )

  extreme_labels <- joined |>
    filter(!is.na(pc1_extreme_label)) |>
    transmute(
      slide_id,
      patient_id,
      PC1,
      label = pc1_extreme_label,
      label_name = pc1_extreme_group
    )

  list(
    median = median_labels,
    extreme = extreme_labels
  )
}

# ------------------------------------------------------------------------------
# Host-gene ranking for GSEA
# ------------------------------------------------------------------------------

run_gene_pc1_ranking <- function(gene_matrix, sample_ids, pc1_values) {
  gene_matrix <- subset_matrix_to_samples(gene_matrix, sample_ids)
  storage.mode(gene_matrix) <- "numeric"

  results <- lapply(
    seq_len(nrow(gene_matrix)),
    function(i) {
      gene_values <- as.numeric(gene_matrix[i, ])
      complete_index <- complete.cases(gene_values, pc1_values)

      if (sum(complete_index) < 10L ||
          length(unique(gene_values[complete_index])) < 2L) {
        return(
          data.frame(
            gene = rownames(gene_matrix)[i],
            rho = NA_real_,
            p_value = NA_real_
          )
        )
      }

      test <- suppressWarnings(
        cor.test(
          gene_values[complete_index],
          pc1_values[complete_index],
          method = "spearman",
          exact = FALSE
        )
      )

      data.frame(
        gene = rownames(gene_matrix)[i],
        rho = unname(test$estimate),
        p_value = test$p.value
      )
    }
  )

  bind_rows(results) |>
    mutate(
      FDR = p.adjust(p_value, method = "BH"),

      # PCA sign is oriented so positive PC1 corresponds to the SINE/Alu side.
      # The thesis GSEA display used the opposite sign so positive enrichment
      # corresponds to the low-PC1 direction.
      rank_low_PC1_positive = -rho
    ) |>
    arrange(desc(rank_low_PC1_positive))
}

# ------------------------------------------------------------------------------
# Optional HoVer-Net integration
# ------------------------------------------------------------------------------

summarize_hovernet_by_pc1_group <- function(
    hovernet_table,
    slide_labels,
    group_column = "pc1_median_group"
) {
  slide_id_col <- first_existing_column(
    hovernet_table,
    c("slide_id", "slide", "filename"),
    "HoVer-Net slide ID"
  )

  joined <- hovernet_table |>
    mutate(slide_id = as.character(.data[[slide_id_col]])) |>
    inner_join(
      slide_labels |>
        select(slide_id, patient_id, label_name),
      by = "slide_id"
    ) |>
    rename(pc1_group = label_name)

  composition_columns <- colnames(joined)[
    grepl("_prop$|immune_epithelial_ratio$", colnames(joined))
  ]

  group_summary <- joined |>
    group_by(pc1_group) |>
    summarise(
      across(
        all_of(composition_columns),
        list(
          median = ~ median(.x, na.rm = TRUE),
          mean = ~ mean(.x, na.rm = TRUE)
        )
      ),
      .groups = "drop"
    )

  list(
    joined = joined,
    summary = group_summary
  )
}

# ==============================================================================
# Main workflow
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load clinical and TE-expression data
# ------------------------------------------------------------------------------

clinical <- readRDS(clinical_path)
te_matrix <- readRDS(te_matrix_path)

if (!is.data.frame(clinical)) {
  stop("The clinical RDS must contain a data frame.")
}

if (!is.matrix(te_matrix) && !is.data.frame(te_matrix)) {
  stop("The TE-expression RDS must contain a matrix or data frame.")
}

te_matrix <- as.matrix(te_matrix)
storage.mode(te_matrix) <- "numeric"

if (is.null(rownames(te_matrix))) {
  stop("The TE-expression matrix must have TE names as row names.")
}

clinical <- align_clinical_to_matrix(clinical, te_matrix)

msi_col <- first_existing_column(
  clinical,
  c("MSI_status", "msi_status", "MSI", "msi"),
  "MSI status"
)

seq_col <- first_existing_column(
  clinical,
  c("seq_type", "sequencing_type"),
  "Sequencing type"
)

cohort_index <- is_mss_value(clinical[[msi_col]]) &
  is_paired_end_value(clinical[[seq_col]])

clinical_mss_pe <- clinical[cohort_index, , drop = FALSE]
te_mss_pe <- te_matrix[, clinical_mss_pe$sample_id, drop = FALSE]

if (nrow(clinical_mss_pe) == 0L) {
  stop("No MSS paired-end samples remained after filtering.")
}

cat("MSS paired-end cohort:", nrow(clinical_mss_pe), "samples\n")
cat("TE features before variance filtering:", nrow(te_mss_pe), "\n")

write.csv(
  clinical_mss_pe |>
    select(sample_id, patient_id, all_of(msi_col), all_of(seq_col)),
  file.path(output_dir, "mss_pe_cohort_manifest.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 2. PCA of TE expression
# ------------------------------------------------------------------------------

pca_result <- run_te_pca(te_mss_pe)

pc1_table <- clinical_mss_pe
pc1_table$PC1 <- pca_result$pca$x[, 1]
pc1_table$PC2 <- pca_result$pca$x[, 2]

variance_table <- data.frame(
  component = paste0("PC", seq_along(pca_result$variance_explained)),
  variance_explained = pca_result$variance_explained,
  cumulative_variance = cumsum(pca_result$variance_explained)
)

loading_summary <- summarize_pc1_loadings(pca_result$pca$rotation)

write.csv(
  pc1_table,
  file.path(output_dir, "mss_pe_pca_scores.csv"),
  row.names = FALSE
)

write.csv(
  variance_table,
  file.path(output_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

write.csv(
  loading_summary$loading_table,
  file.path(output_dir, "pc1_te_loadings.csv"),
  row.names = FALSE
)

write.csv(
  loading_summary$family_summary,
  file.path(output_dir, "pc1_te_family_summary.csv"),
  row.names = FALSE
)

writeLines(
  paste("PC1 sign flipped for interpretation:", pca_result$pc1_flipped),
  file.path(output_dir, "pc1_orientation_note.txt")
)

# ------------------------------------------------------------------------------
# 3. Survival analysis
# ------------------------------------------------------------------------------

survival_result <- run_pc1_survival(pc1_table)

write.csv(
  survival_result$continuous,
  file.path(output_dir, "pc1_continuous_cox_model.csv"),
  row.names = FALSE
)

write.csv(
  survival_result$logrank,
  file.path(output_dir, "pc1_median_logrank.csv"),
  row.names = FALSE
)

# The original exploratory analysis also applied KAPS to PC1. That step can use
# the KAPS workflow documented in Part 1 once the historical package environment
# has been verified.

# ------------------------------------------------------------------------------
# 4. Generate PC1 groups for downstream WSI analysis
# ------------------------------------------------------------------------------

pc1_groups <- make_pc1_groups(
  pc1_table |>
    select(sample_id, patient_id, PC1)
)

write.csv(
  pc1_groups,
  file.path(output_dir, "pc1_patient_groups.csv"),
  row.names = FALSE
)

if (path_is_available(slide_metadata_path)) {
  slide_metadata <- read.csv(
    slide_metadata_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  clam_labels <- prepare_slide_labels(pc1_groups, slide_metadata)

  write.csv(
    clam_labels$median,
    file.path(output_dir, "clam_labels_pc1_median.csv"),
    row.names = FALSE
  )

  write.csv(
    clam_labels$extreme,
    file.path(output_dir, "clam_labels_pc1_extreme_quartiles.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Immune and stromal correlations
# ------------------------------------------------------------------------------

if (path_is_available(immune_scores_path)) {
  immune_scores <- readRDS(immune_scores_path)
  aligned_scores <- align_sample_table(
    immune_scores,
    pc1_table$sample_id
  )

  immune_correlations <- run_pc1_score_correlations(
    pc1_table = pc1_table,
    score_table = aligned_scores
  )

  write.csv(
    immune_correlations,
    file.path(output_dir, "pc1_immune_stromal_correlations.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 6. Host-gene ranking for GSEA
# ------------------------------------------------------------------------------

if (path_is_available(gene_expression_path)) {
  gene_matrix <- readRDS(gene_expression_path)

  if (!is.matrix(gene_matrix) && !is.data.frame(gene_matrix)) {
    stop("The host-gene expression RDS must contain a matrix or data frame.")
  }

  gene_matrix <- as.matrix(gene_matrix)

  if (is.null(rownames(gene_matrix))) {
    stop("The gene-expression matrix must have gene names as row names.")
  }

  gene_ranking <- run_gene_pc1_ranking(
    gene_matrix = gene_matrix,
    sample_ids = pc1_table$sample_id,
    pc1_values = pc1_table$PC1
  )

  write.csv(
    gene_ranking,
    file.path(output_dir, "gene_pc1_correlations.csv"),
    row.names = FALSE
  )

  # Two-column preranked GSEA input: gene and ranking statistic.
  write.table(
    gene_ranking |>
      filter(!is.na(rank_low_PC1_positive)) |>
      select(gene, rank_low_PC1_positive),
    file.path(output_dir, "gene_pc1_low_direction.rnk"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 7. Optional HoVer-Net cell-composition comparison
# ------------------------------------------------------------------------------

if (
  path_is_available(hovernet_composition_path) &&
  exists("clam_labels") &&
  nrow(clam_labels$median) > 0L
) {
  hovernet_table <- read.csv(
    hovernet_composition_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  hovernet_result <- summarize_hovernet_by_pc1_group(
    hovernet_table = hovernet_table,
    slide_labels = clam_labels$median
  )

  write.csv(
    hovernet_result$joined,
    file.path(output_dir, "hovernet_pc1_slide_composition.csv"),
    row.names = FALSE
  )

  write.csv(
    hovernet_result$summary,
    file.path(output_dir, "hovernet_pc1_group_summary.csv"),
    row.names = FALSE
  )
}

cat("\nPart 3 workflow complete.\n")
cat("Outputs written to:", normalizePath(output_dir), "\n")
