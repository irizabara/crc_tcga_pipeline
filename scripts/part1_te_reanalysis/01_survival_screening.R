#!/usr/bin/env Rscript
# Cox Regression Analysis: Test Each Individual TE
# Three model types for each endpoint (OS, DSS, DFI, PFI)

library(survival)
library(dplyr)
library(tidyr)
library(progress)

# Load data
clinical <- readRDS("~/nas/irizab/TE_clusters/02_clusters/clinical_final.RDS")
te_matrix <- readRDS("~/nas/irizab/TE_clusters/02_clusters/te_matrix_matched_614samples.RDS")

cat("=================================================================\n")
cat("STEP 1: DATA PREPARATION\n")
cat("=================================================================\n\n")

# ============================================================================
# Prepare clinical data with proper variable types
# ============================================================================

# Convert age to numeric
clinical$age <- as.numeric(clinical$age_at_initial_pathologic_diagnosis)

# Clean and simplify stage
clinical$stage_raw <- clinical$ajcc_pathologic_tumor_stage
clinical$stage <- case_when(
  grepl("Stage I[^I]|Stage I$", clinical$stage_raw) ~ "I",
  grepl("Stage II", clinical$stage_raw) ~ "II",
  grepl("Stage III", clinical$stage_raw) ~ "III",
  grepl("Stage IV", clinical$stage_raw) ~ "IV",
  TRUE ~ NA_character_
)
clinical$stage <- factor(clinical$stage, levels = c("I", "II", "III", "IV"))

# Convert seq_type to factor
clinical$seq_type <- factor(clinical$seq_type)

# Check the variables
cat("Clinical variable summary:\n")
cat("Age range:", range(clinical$age, na.rm = TRUE), "\n")
cat("Stage distribution:\n")
print(table(clinical$stage, useNA = "ifany"))
cat("\nSeq_type distribution:\n")
print(table(clinical$seq_type, useNA = "ifany"))
cat("\n")

# Verify sample matching
stopifnot(all(colnames(te_matrix) == paste0(clinical$bcr_patient_barcode, "-01A")))

cat("Data dimensions:\n")
cat("  Clinical samples:", nrow(clinical), "\n")
cat("  TE features:", nrow(te_matrix), "\n")
cat("  Samples in TE matrix:", ncol(te_matrix), "\n\n")

# ============================================================================
# STEP 2: FUNCTION TO RUN THREE MODELS FOR ONE TE
# ============================================================================

run_te_models <- function(te_name, te_expr, clinical_data, endpoint = "OS") {

  # Prepare survival variables
  time_var <- paste0(endpoint, ".time")
  status_var <- endpoint

  # Create analysis dataframe
  analysis_df <- data.frame(
    time = clinical_data[[time_var]],
    status = clinical_data[[status_var]],
    te_expr = te_expr,
    age = clinical_data$age,
    stage = clinical_data$stage,
    seq_type = clinical_data$seq_type
  )

  # Remove rows with missing data
  analysis_df_complete <- analysis_df[complete.cases(analysis_df), ]

  # Check if we have enough data
  if (nrow(analysis_df_complete) < 50 || sum(analysis_df_complete$status) < 10) {
    return(NULL)
  }

  # Dichotomize TE expression at median
  analysis_df_complete$TE_g <- ifelse(
    analysis_df_complete$te_expr > median(analysis_df_complete$te_expr),
    "High",
    "Low"
  )
  analysis_df_complete$TE_g <- factor(analysis_df_complete$TE_g, levels = c("Low", "High"))

  # Initialize results list
  results <- list()

  # ========================================================================
  # MODEL 1: Primary model (biology-adjusted)
  # OS ~ TE_g + age + stage
  # ========================================================================
  tryCatch({
    model1 <- coxph(Surv(time, status) ~ TE_g + age + stage,
                    data = analysis_df_complete)

    coef1 <- summary(model1)$coefficients
    conf1 <- summary(model1)$conf.int

    # Extract TE_g results
    results$model1 <- data.frame(
      TE = te_name,
      endpoint = endpoint,
      model = "Primary (biology-adjusted)",
      formula = "~ TE_g + age + stage",
      n = nrow(analysis_df_complete),
      n_events = sum(analysis_df_complete$status),
      coef = coef1["TE_gHigh", "coef"],
      se = coef1["TE_gHigh", "se(coef)"],
      HR = conf1["TE_gHigh", "exp(coef)"],
      lower_CI = conf1["TE_gHigh", "lower .95"],
      upper_CI = conf1["TE_gHigh", "upper .95"],
      z = coef1["TE_gHigh", "z"],
      p_value = coef1["TE_gHigh", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results$model1 <<- data.frame(
      TE = te_name, endpoint = endpoint,
      model = "Primary (biology-adjusted)",
      formula = "~ TE_g + age + stage",
      n = NA, n_events = NA, coef = NA, se = NA,
      HR = NA, lower_CI = NA, upper_CI = NA, z = NA, p_value = NA,
      error = as.character(e)
    )
  })

  # ========================================================================
  # MODEL 2: Batch sensitivity with stratification (Iriza's favorite!)
  # OS ~ TE_g + age + stage + strata(seq_type)
  # ========================================================================
  tryCatch({
    model2 <- coxph(Surv(time, status) ~ TE_g + age + stage + strata(seq_type),
                    data = analysis_df_complete)

    coef2 <- summary(model2)$coefficients
    conf2 <- summary(model2)$conf.int

    results$model2 <- data.frame(
      TE = te_name,
      endpoint = endpoint,
      model = "Batch sensitivity (stratified)",
      formula = "~ TE_g + age + stage + strata(seq_type)",
      n = nrow(analysis_df_complete),
      n_events = sum(analysis_df_complete$status),
      coef = coef2["TE_gHigh", "coef"],
      se = coef2["TE_gHigh", "se(coef)"],
      HR = conf2["TE_gHigh", "exp(coef)"],
      lower_CI = conf2["TE_gHigh", "lower .95"],
      upper_CI = conf2["TE_gHigh", "upper .95"],
      z = coef2["TE_gHigh", "z"],
      p_value = coef2["TE_gHigh", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results$model2 <<- data.frame(
      TE = te_name, endpoint = endpoint,
      model = "Batch sensitivity (stratified)",
      formula = "~ TE_g + age + stage + strata(seq_type)",
      n = NA, n_events = NA, coef = NA, se = NA,
      HR = NA, lower_CI = NA, upper_CI = NA, z = NA, p_value = NA,
      error = as.character(e)
    )
  })

  # ========================================================================
  # MODEL 3: Additional sensitivity with seq_type as covariate
  # OS ~ TE_g + age + stage + seq_type
  # ========================================================================
  tryCatch({
    model3 <- coxph(Surv(time, status) ~ TE_g + age + stage + seq_type,
                    data = analysis_df_complete)

    coef3 <- summary(model3)$coefficients
    conf3 <- summary(model3)$conf.int

    results$model3 <- data.frame(
      TE = te_name,
      endpoint = endpoint,
      model = "Additional sensitivity (covariate)",
      formula = "~ TE_g + age + stage + seq_type",
      n = nrow(analysis_df_complete),
      n_events = sum(analysis_df_complete$status),
      coef = coef3["TE_gHigh", "coef"],
      se = coef3["TE_gHigh", "se(coef)"],
      HR = conf3["TE_gHigh", "exp(coef)"],
      lower_CI = conf3["TE_gHigh", "lower .95"],
      upper_CI = conf3["TE_gHigh", "upper .95"],
      z = coef3["TE_gHigh", "z"],
      p_value = coef3["TE_gHigh", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results$model3 <<- data.frame(
      TE = te_name, endpoint = endpoint,
      model = "Additional sensitivity (covariate)",
      formula = "~ TE_g + age + stage + seq_type",
      n = NA, n_events = NA, coef = NA, se = NA,
      HR = NA, lower_CI = NA, upper_CI = NA, z = NA, p_value = NA,
      error = as.character(e)
    )
  })

  # Combine all three models
  all_results <- bind_rows(results)
  return(all_results)
}

# ============================================================================
# STEP 3: RUN MODELS FOR ALL TEs AND ALL ENDPOINTS
# ============================================================================

cat("=================================================================\n")
cat("STEP 2: RUNNING COX MODELS\n")
cat("=================================================================\n\n")

cat("Configuration:\n")
cat("  - Number of TEs: ", nrow(te_matrix), "\n")
cat("  - Number of endpoints: 4 (OS, DSS, DFI, PFI)\n")
cat("  - Models per TE: 3 (Primary, Stratified, Covariate)\n")
cat("  - Total models to run: ", nrow(te_matrix) * 4 * 3, "\n\n")

cat("Time estimate:\n")
cat("  - If 1 model took 2 minutes for all TEs\n")
cat("  - Then 3 models will take ~6 minutes per endpoint\n")
cat("  - Total estimated time: ~24 minutes for all 4 endpoints\n\n")
cat("Starting analysis...\n")

overall_start <- Sys.time()

endpoints <- c("OS", "DSS", "DFI", "PFI")
all_results_list <- list()

for (endpoint in endpoints) {
  cat(paste0("\n>>> Processing endpoint: ", endpoint, " <<<\n"))
  cat(paste0("Total TEs to process: ", nrow(te_matrix), "\n"))
  cat("Running 3 models per TE (Primary, Stratified, Covariate)\n\n")

  # Create progress bar
  pb <- progress_bar$new(
    format = "  [:bar] :percent | :current/:total TEs | Elapsed: :elapsed | ETA: :eta",
    total = nrow(te_matrix),
    clear = FALSE,
    width = 80
  )

  endpoint_results <- list()
  start_time <- Sys.time()

  for (i in 1:nrow(te_matrix)) {

    # Update progress bar
    pb$tick()

    te_name <- rownames(te_matrix)[i]
    te_expr <- as.numeric(te_matrix[i, ])

    # Run models for this TE
    te_results <- run_te_models(
      te_name = te_name,
      te_expr = te_expr,
      clinical_data = clinical,
      endpoint = endpoint
    )

    if (!is.null(te_results)) {
      endpoint_results[[i]] <- te_results
    }
  }

  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "mins"))

  # Combine results for this endpoint
  endpoint_df <- bind_rows(endpoint_results)
  all_results_list[[endpoint]] <- endpoint_df

  cat(paste0("\n  ✓ Completed ", endpoint, ": ", nrow(endpoint_df), " results in ",
             round(time_taken, 2), " minutes\n"))
}

# Combine all results
cat("\n=================================================================\n")
cat("STEP 3: COMBINING AND SAVING RESULTS\n")
cat("=================================================================\n\n")

overall_end <- Sys.time()
total_time <- as.numeric(difftime(overall_end, overall_start, units = "mins"))
cat("Total analysis time:", round(total_time, 2), "minutes\n\n")

final_results <- bind_rows(all_results_list)

# Add FDR correction for p-values (within each model type and endpoint)
final_results <- final_results %>%
  group_by(endpoint, model) %>%
  mutate(
    fdr = p.adjust(p_value, method = "fdr"),
    bonferroni = p.adjust(p_value, method = "bonferroni")
  ) %>%
  ungroup()

# Summary statistics
cat("Total results generated:", nrow(final_results), "\n")
cat("\nBreakdown by endpoint and model:\n")
print(table(final_results$endpoint, final_results$model))

cat("\nSignificant results (p < 0.05):\n")
sig_summary <- final_results %>%
  filter(!is.na(p_value)) %>%
  group_by(endpoint, model) %>%
  summarise(
    total = n(),
    sig_p05 = sum(p_value < 0.05),
    sig_fdr05 = sum(fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  )
print(sig_summary)

# Save results
output_file <- "~/nas/irizab/TE_clusters/02_clusters/cox_results_all_TEs_three_models.RDS"
saveRDS(final_results, output_file)
cat("\nResults saved to:", output_file, "\n")

# Also save as CSV for easy viewing
csv_file <- "~/nas/irizab/TE_clusters/02_clusters/cox_results_all_TEs_three_models.csv"
write.csv(final_results, csv_file, row.names = FALSE)
cat("Results also saved as CSV:", csv_file, "\n")

# ============================================================================
# STEP 4: CREATE SUMMARY TABLES
# ============================================================================

cat("\n=================================================================\n")
cat("STEP 4: CREATING SUMMARY TABLES\n")
cat("=================================================================\n\n")

# Top significant TEs for each model and endpoint
top_tes <- final_results %>%
  filter(!is.na(p_value)) %>%
  group_by(endpoint, model) %>%
  arrange(p_value) %>%
  slice_head(n = 20) %>%
  ungroup()

top_file <- "~/nas/irizab/TE_clusters/02_clusters/cox_top20_TEs_per_model.csv"
write.csv(top_tes, top_file, row.names = FALSE)
cat("Top 20 TEs per model saved to:", top_file, "\n")

# Print top results for OS with primary model
cat("\n--- Top 10 TEs for OS (Primary model) ---\n")
os_primary_top <- final_results %>%
  filter(endpoint == "OS", model == "Primary (biology-adjusted)", !is.na(p_value)) %>%
  arrange(p_value) %>%
  select(TE, HR, lower_CI, upper_CI, p_value, fdr) %>%
  head(10)
print(os_primary_top)

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=================================================================\n")
cat("\nFiles generated:\n")
cat("1. Full results (RDS):", output_file, "\n")
cat("2. Full results (CSV):", csv_file, "\n")
cat("3. Top 20 TEs per model:", top_file, "\n")