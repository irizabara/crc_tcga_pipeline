# Confounder-Aware Analysis of Transposable Element Expression, Immune Heterogeneity, and Whole-Slide Morphology in Colorectal Cancer

This repository accompanies the MSc thesis:

> **Confounder-Aware Analysis of Transposable Element Expression, Immune Heterogeneity, and Whole-Slide Morphology in Colorectal Cancer**  
> **Iriza Baranyanka**  
> MSc Bioinformatics and Systems Biology, Vrije Universiteit Amsterdam and University of Amsterdam, 2026

## Overview

This project investigated whether transposable element (TE) expression states in colorectal cancer are:

1. robust to technical and clinical confounding,
2. associated with survival and immune heterogeneity, and
3. reflected in H&E whole-slide morphology.

The analysis used TCGA colorectal cancer samples from the COAD and READ cohorts. The matched molecular cohort contained 614 primary tumour samples, including 351 paired-end and 263 single-end RNA-seq samples.

This repository is **not a complete archive of every exploratory script produced during the project**. It is a curated thesis companion containing representative, readable workflows that show what was done, why each step was performed, what inputs were expected, and how the main outputs were generated.

## Study workflow

![Overview of the thesis analysis workflow](docs/workflow_overview.png)

## Repository structure

```text
crc_tcga_pipeline/
├── README.md
├── docs/
│   ├── environment_r.yml
│   ├── environment_kaps.yml
│   ├── kaps_package_description.txt
│   └── workflow_overview.png
└── scripts/
    ├── part1_te_reanalysis/
    │   ├── 01_survival_screening.R
    │   ├── 02_immune_correlation_screening.R
    │   └── 03_te_score_and_kaps.R
    ├── part2_wsi_prediction/
    │   ├── 01_clam_wsi_pipeline.py
    │   ├── 02_hovernet_cell_composition.py
    │   └── part2_clam_hovernet_workflow.md
    └── part3_mss_pc1/
        └── 01_mss_pe_pc1_analysis.R
```

## Part 1 — Confounder-aware TE reanalysis

Part 1 evaluated whether TE associations with outcome and immune variation remained interpretable after accounting for technical and clinical sources of variation.

### 1. Survival screening

[`01_survival_screening.R`](scripts/part1_te_reanalysis/01_survival_screening.R)

This workflow performs TE-by-TE Cox regression across four endpoints:

- overall survival,
- disease-specific survival,
- disease-free interval, and
- progression-free interval.

For each TE, expression is dichotomised at the cohort median. The script compares:

- a primary model adjusted for age and pathological stage, and
- a sequencing-sensitive model stratified by sequencing type.

The final selection rule retains TEs with support in at least two primary-model endpoints and excludes TEs showing a hazard-ratio direction flip between the primary and sequencing-stratified models in any evaluable endpoint.

The purpose of the sequencing analysis is not to reproduce a sequencing effect as a biological finding, but to demonstrate how study-specific technical variation was tested and controlled.

### 2. TE–immune correlation screening

[`02_immune_correlation_screening.R`](scripts/part1_te_reanalysis/02_immune_correlation_screening.R)

This workflow evaluates associations between TE expression and immune indices after confounder adjustment.

Two analyses are represented:

- primary residualisation using age and pathological stage;
- sensitivity residualisation using age, pathological stage, and sequencing type.

Spearman correlations are calculated between residualised TE expression and residualised immune scores. Benjamini-Hochberg correction is applied across the immune-index panel for each TE.

The script produces:

- full primary and sensitivity correlation tables,
- per-TE immune association summaries,
- direction-concordance checks, and
- a set of robust immune-associated candidate TEs.

### 3. TE-score generation and KAPS

[`03_te_score_and_kaps.R`](scripts/part1_te_reanalysis/03_te_score_and_kaps.R)

This workflow documents TE-score generation and K-adaptive partitioning for survival analysis.

The final seven-TE signature was:

- `LTR106_Mam`
- `LTR19-int`
- `LTR80A`
- `MER44B`
- `MER57E1`
- `MER65C`
- `Tigger11a`

The score was generated from the normalized expression of the selected TEs and evaluated using KAPS-based survival partitioning. The script also documents the comparison with the published Zhu nine-TE signature.

The original KAPS analysis was performed in an institutional server environment. The current script was reconstructed from the final thesis methods, retained outputs, and project history. It is included to document the analytical logic rather than to claim byte-for-byte identity with the historical server-side script.

## Part 2 — Whole-slide image prediction with CLAM and HoVer-Net

Part 2 examined whether TE-derived molecular states could be predicted from H&E whole-slide images.

### 1. CLAM workflow skeleton

[`01_clam_wsi_pipeline.py`](scripts/part2_wsi_prediction/01_clam_wsi_pipeline.py)

This curated template documents the main CLAM workflow:

1. prepare slide-level labels,
2. preserve patient-level separation between data splits,
3. segment tissue and extract patches,
4. generate patch embeddings using a pretrained UNI2-h encoder,
5. train and evaluate CLAM,
6. combine fold-level predictions.

The script is intentionally a readable workflow skeleton rather than a replacement for the external CLAM repository. Local paths, the exact CLAM fork, model weights, and source slides must be supplied by the user.

### 2. HoVer-Net cell-composition summary

[`02_hovernet_cell_composition.py`](scripts/part2_wsi_prediction/02_hovernet_cell_composition.py)

This workflow represents the downstream analysis performed after HoVer-Net inference.

It:

- reads per-slide HoVer-Net nuclei predictions,
- counts epithelial, lymphocyte, macrophage, neutrophil, and unlabeled nuclei,
- calculates cell-type proportions,
- derives a combined immune-cell proportion,
- calculates an immune-to-epithelial ratio, and
- optionally joins cell-composition summaries to study groups such as low- and high-PC1 labels.

### 3. CLAM and HoVer-Net workflow notes

[`part2_clam_hovernet_workflow.md`](scripts/part2_wsi_prediction/part2_clam_hovernet_workflow.md)

This tutorial documents the external CLAM, UNI2-h, and HoVer-Net workflow used during the project. It was adapted from an internal guide prepared by **Mark Yang** and is shared with permission.

The document includes:

- whole-slide patch extraction,
- UNI2-h feature extraction,
- HoVer-Net inference,
- nuclei-class definitions, and
- optional integration of CLAM attention with cell-composition summaries.

Institution-specific paths, credentials, and access tokens were removed.

## Part 3 — MSS paired-end PC1 analysis

[`01_mss_pe_pc1_analysis.R`](scripts/part3_mss_pc1/01_mss_pe_pc1_analysis.R)

Part 3 narrowed the analysis to microsatellite-stable, paired-end samples to reduce heterogeneity and examine a sequencing-qualified TE-expression axis.

The workflow documents:

1. filtering to MSS paired-end samples,
2. aligning clinical and TE-expression data,
3. performing PCA on normalized TE expression,
4. orienting and interpreting PC1 using TE-family loadings,
5. summarizing SINE/Alu, LTR/ERV, LINE/L1/L2, MER, and DNA-transposon contributions,
6. testing PC1 associations with survival,
7. testing PC1 correlations with immune and stromal indices,
8. generating median-split and extreme-quartile labels for CLAM,
9. producing a host-gene ranking for GSEA, and
10. optionally integrating HoVer-Net cell-composition summaries.

This part treats PC1 as a TE-family composition axis rather than a simple measure of total TE burden.

## Inputs and outputs

The scripts expect locally prepared analysis objects rather than raw TCGA downloads.

Typical inputs include:

- normalized TE-expression matrices,
- clinical and survival annotations,
- immune or stromal score tables,
- host-gene expression matrices,
- slide-level metadata,
- CLAM feature tensors,
- CLAM prediction outputs, and
- HoVer-Net JSON outputs.

Typical outputs include:

- Cox regression tables,
- survival candidate lists,
- TE–immune correlation tables,
- TE-score and KAPS summaries,
- PCA scores and loadings,
- PC1-derived CLAM label tables,
- ranked gene lists for GSEA, and
- per-slide cell-composition summaries.

## Data availability

This repository does **not** redistribute:

- TCGA RNA-seq data,
- patient-level clinical or molecular metadata,
- whole-slide images,
- extracted patch embeddings,
- model checkpoints,
- patient-level predictions,
- large intermediate RDS objects, or
- institutional credentials and filesystem paths.

Users must obtain TCGA data through the appropriate data portals and access procedures, then adapt the scripts to their local file organization.

## External tools and resources

The project used or referenced:

- REdiscoverTE
- Salmon
- CLAM
- UNI / UNI2-h
- HoVer-Net
- KAPS
- MSigDB Hallmark gene sets
- TCGA COAD and READ data resources

The REdiscoverTE implementation used during the project is not redistributed because the associated public repositories and installation links are no longer available.

External source code, pretrained weights, and datasets remain subject to their original licences and access conditions.

## Environments

Two environment specifications are included:

- [`environment_r.yml`](docs/environment_r.yml) — general R analysis environment
- [`environment_kaps.yml`](docs/environment_kaps.yml) — KAPS-specific environment

The exact KAPS package metadata is recorded in:

- [`kaps_package_description.txt`](docs/kaps_package_description.txt)

The recovered KAPS installation corresponds to package version `1.1.4`, built from GitHub ref `v1.1.5`.

Environment files are historical snapshots and may require minor channel or version adjustments on a new system.

To recreate an environment:

```bash
conda env create -f docs/environment_r.yml
```

or:

```bash
conda env create -f docs/environment_kaps.yml
```

## How to use this repository

Clone the repository:

```bash
git clone git@github.com:irizabara/crc_tcga_pipeline.git
cd crc_tcga_pipeline
```

Then read the workflow in order:

1. `scripts/part1_te_reanalysis/`
2. `scripts/part2_wsi_prediction/`
3. `scripts/part3_mss_pc1/`

The scripts are designed to be read and adapted. They are not guaranteed to run unchanged because the original project depended on institutional compute infrastructure, non-redistributable data, external repositories, and pretrained models.

Where command-line arguments are supported, local input and output paths can be supplied directly. Other scripts contain clearly marked placeholder paths that should be replaced before use.

## Reproducibility and FAIR principles

The repository makes the major analytical decisions inspectable:

- cohort restriction,
- confounder handling,
- model formulas,
- survival endpoints,
- TE-selection criteria,
- immune-correlation logic,
- TE-score construction,
- WSI label preparation,
- PC1 interpretation, and
- downstream image and cell-composition analysis.

## Important limitations

- Some original analyses were run on institutional servers that are no longer accessible.
- Some scripts were cleaned or reconstructed after the analyses were completed.
- Exact numerical reproduction may depend on package versions, preprocessing details, random seeds, local code modifications, and the availability of external model weights.

**Iriza Baranyanka**  
GitHub: [@irizabara](https://github.com/irizabara)
