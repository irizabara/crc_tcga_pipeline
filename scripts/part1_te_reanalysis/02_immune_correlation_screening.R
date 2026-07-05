# ==============================================================================
# TE-Immune Correlation Analysis with Proper Confounding Control
# ==============================================================================

library(dplyr)
library(tidyr)

# ==============================================================================
# STEP 0: Sanity Checks
# ==============================================================================

cat("=== Step 0: Sanity Checks ===\n")

# Check dimensions
stopifnot(nrow(immune) == nrow(clinical))
stopifnot(ncol(te_matrix) == nrow(clinical))

# Check sample ID alignment
stopifnot(all(colnames(te_matrix) == paste0(clinical$bcr_patient_barcode, "-01A")))

cat("✓ All samples aligned: n =", nrow(clinical), "\n")
cat("✓ TEs:", nrow(te_matrix), "\n")
cat("✓ Immune indices:", ncol(immune), "\n\n")

# ==============================================================================
# STEP 1: Residualization Helper Function
# ==============================================================================

residualize <- function(y, covariates_df) {
  # Build complete data frame
  df <- data.frame(y = y, covariates_df)

  # Identify complete cases
  complete_idx <- complete.cases(df)

  # Initialize output (NAs where input was incomplete)
  resid_out <- rep(NA_real_, length(y))

  # Fit model on complete cases only
  if (sum(complete_idx) > 0) {
    fit <- lm(y ~ ., data = df[complete_idx, , drop = FALSE])
    resid_out[complete_idx] <- residuals(fit)
  }

  return(resid_out)
}

# ==============================================================================
# STEP 2 & 3: Residualization Wrapper
# ==============================================================================

run_residualization <- function(te_matrix, immune, clinical,
                                include_seq_type = FALSE) {

  cat("=== Residualizing with covariates:",
      ifelse(include_seq_type, "age + stage + seq_type", "age + stage"),
      "===\n")

  # Define covariate data frame
  if (include_seq_type) {
    covariates <- clinical[, c("age", "stage", "seq_type")]
  } else {
    covariates <- clinical[, c("age", "stage")]
  }

  # STEP 2: Residualize immune (29 columns)
  cat("Residualizing immune indices...\n")
  immune_resid <- immune
  for (j in 1:ncol(immune)) {
    immune_resid[, j] <- residualize(immune[, j], covariates)
  }

  # STEP 3: Residualize TEs (1076 rows)
  cat("Residualizing", nrow(te_matrix), "TEs...\n")
  te_resid <- te_matrix  # Initialize with same structure
  for (i in 1:nrow(te_matrix)) {
    te_vec <- as.numeric(te_matrix[i, ])  # Extract TE across samples
    te_resid[i, ] <- residualize(te_vec, covariates)
  }

  cat("✓ Residualization complete\n\n")

  return(list(
    te_resid = te_resid,
    immune_resid = immune_resid
  ))
}

# ==============================================================================
# STEP 4 & 5: Correlation + Per-TE FDR
# ==============================================================================

compute_correlations <- function(te_resid, immune_resid, method = "spearman") {

  cat("=== Step 4: Computing correlations ===\n")

  # Initialize results list
  results_list <- list()

  # Loop through TEs
  for (i in 1:nrow(te_resid)) {
    te_name <- rownames(te_resid)[i]
    te_vec <- as.numeric(te_resid[i, ])

    # Store correlations for this TE across all immune indices
    te_results <- data.frame(
      TE = character(ncol(immune_resid)),
      immune_index = character(ncol(immune_resid)),
      n = integer(ncol(immune_resid)),
      rho = numeric(ncol(immune_resid)),
      p_value = numeric(ncol(immune_resid)),
      stringsAsFactors = FALSE
    )

    # Loop through immune indices
    for (j in 1:ncol(immune_resid)) {
      immune_name <- colnames(immune_resid)[j]
      immune_vec <- immune_resid[, j]

      # Correlation test (handles NAs automatically with complete.obs)
      test <- cor.test(te_vec, immune_vec,
                      method = method,
                      use = "complete.obs",
                      exact = FALSE)  # For ties in Spearman

      # Store results
      te_results$TE[j] <- te_name
      te_results$immune_index[j] <- immune_name
      te_results$n[j] <- sum(complete.cases(te_vec, immune_vec))
      te_results$rho[j] <- test$estimate
      te_results$p_value[j] <- test$p.value
    }

    # STEP 5: Per-TE FDR correction (across 29 immune indices)
    te_results$FDR <- p.adjust(te_results$p_value, method = "BH")

    results_list[[i]] <- te_results

    # Progress
    if (i %% 100 == 0) cat("  Processed", i, "TEs...\n")
  }

  # Combine all results
  full_table <- bind_rows(results_list)

  cat("✓ Correlation analysis complete:", nrow(full_table), "tests\n\n")

  return(full_table)
}

# ==============================================================================
# STEP 6: Define Immunogenic TEs
# ==============================================================================

summarize_immunogenic_TEs <- function(full_table,
                                     rho_threshold = 0.4,
                                     p_threshold = 1e-4,
                                     fdr_threshold = 0.05) {

  cat("=== Step 6: Identifying Immunogenic TEs ===\n")

  # Apply thresholds to identify significant associations
  full_table <- full_table %>%
    mutate(
      sig_zhu = abs(rho) >= rho_threshold & p_value < p_threshold,
      sig_fdr = FDR < fdr_threshold
    )

  # Summarize per TE
  te_summary <- full_table %>%
    group_by(TE) %>%
    summarise(
      max_abs_rho = max(abs(rho), na.rm = TRUE),
      min_p = min(p_value, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE),
      n_sig_zhu = sum(sig_zhu, na.rm = TRUE),
      n_sig_fdr = sum(sig_fdr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      immunogenic_zhu = max_abs_rho >= rho_threshold & min_p < p_threshold,
      immunogenic_fdr = n_sig_fdr > 0
    ) %>%
    arrange(desc(max_abs_rho))

  # Report counts
  cat("Immunogenic TEs (Zhu criteria: |rho|≥0.4, p<1e-4):",
      sum(te_summary$immunogenic_zhu), "\n")
  cat("Immunogenic TEs (FDR<0.05 for ≥1 immune index):",
      sum(te_summary$immunogenic_fdr), "\n\n")

  return(list(
    full_table = full_table,
    te_summary = te_summary
  ))
}

# ==============================================================================
# MAIN ANALYSIS WRAPPER
# ==============================================================================

run_full_analysis <- function(te_matrix, immune, clinical,
                              include_seq_type = FALSE,
                              label = "") {

  cat("\n")
  cat("################################################################################\n")
  cat("# ANALYSIS:", label, "\n")
  cat("################################################################################\n\n")

  # Residualize
  resid_data <- run_residualization(te_matrix, immune, clinical, include_seq_type)

  # Correlations
  full_table <- compute_correlations(resid_data$te_resid, resid_data$immune_resid)

  # Summarize
  summary_results <- summarize_immunogenic_TEs(full_table)

  return(summary_results)
}

# ==============================================================================
# RUN BOTH ANALYSES
# ==============================================================================

# Analysis 1: age + stage (PRIMARY)
results_primary <- run_full_analysis(
  te_matrix, immune, clinical,
  include_seq_type = FALSE,
  label = "PRIMARY (age + stage)"
)

# Analysis 2: age + stage + seq_type (SENSITIVITY)
results_sensitivity <- run_full_analysis(
  te_matrix, immune, clinical,
  include_seq_type = TRUE,
  label = "SENSITIVITY (age + stage + seq_type)"
)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("=== Saving results ===\n")

# Save full tables
write.csv(results_primary$full_table,
          "te_immune_correlation_full_primary.csv",
          row.names = FALSE)
write.csv(results_sensitivity$full_table,
          "te_immune_correlation_full_sensitivity.csv",
          row.names = FALSE)

# Save TE summaries
write.csv(results_primary$te_summary,
          "te_immune_summary_primary.csv",
          row.names = FALSE)
write.csv(results_sensitivity$te_summary,
          "te_immune_summary_sensitivity.csv",
          row.names = FALSE)

cat("✓ Results saved\n\n")

# ==============================================================================
# QUICK EXPLORATION
# ==============================================================================

cat("=== TOP 20 IMMUNOGENIC TEs (PRIMARY) ===\n")
print(head(results_primary$te_summary %>%
             filter(immunogenic_fdr) %>%
             select(TE, max_abs_rho, min_FDR, n_sig_fdr),
           20))

cat("\n=== OVERLAP CHECK ===\n")
primary_TEs <- results_primary$te_summary %>%
  filter(immunogenic_fdr) %>%
  pull(TE)

sensitivity_TEs <- results_sensitivity$te_summary %>%
  filter(immunogenic_fdr) %>%
  pull(TE)

cat("Primary immunogenic TEs:", length(primary_TEs), "\n")
cat("Sensitivity immunogenic TEs:", length(sensitivity_TEs), "\n")
cat("Overlap:", length(intersect(primary_TEs, sensitivity_TEs)), "\n")