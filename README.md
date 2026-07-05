# Confounder-Aware Analysis of Transposable Element Expression in Colorectal Cancer

Code supporting the MSc thesis:

**Confounder-Aware Analysis of Transposable Element Expression, Immune Heterogeneity, and Whole-Slide Morphology in Colorectal Cancer**

This repository contains analysis scripts used to investigate whether transposable-element expression states in colorectal cancer are robust, biologically interpretable, and recoverable from H&E whole-slide morphology.

## Repository structure

- `scripts/part1_te_reanalysis/`  
  Confounder-aware TE preprocessing, survival screening, immune-correlation analysis, signature derivation, and KAPS analysis.

- `scripts/part2_wsi_prediction/`  
  Label preparation, patient-level cross-validation, CLAM evaluation, MSI positive-control analysis, attention extraction, and HoVer-Net summaries.

- `scripts/part3_mss_pc1/`  
  MSS paired-end PCA, TE-family interpretation, immune and stromal correlations, survival analysis, WSI label preparation, and GSEA ranking.

- `metadata/`  
  Documentation of expected input formats and identifier matching.

- `docs/`  
  Additional workflow and implementation notes.

## Data availability

TCGA RNA-seq data, clinical annotations, whole-slide images, model checkpoints, extracted feature tensors, and patient-level outputs are not redistributed in this repository.

Access to controlled or large TCGA data must be obtained through the appropriate data-access channels.

## External tools

The analyses rely on external tools and resources including:

- REdiscoverTE
- CLAM
- UNI / UNI2-h
- HoVer-Net
- KAPS
- MSigDB Hallmark gene sets

External software is referenced but not redistributed here.

## Reproducibility

Many scripts were developed within an institutional computing environment and contain paths that must be adapted before use. The repository is intended to document the final analytical workflow and the project-specific scripts used in the thesis.
