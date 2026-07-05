# CLAM and HoVer-Net Workflow Notes

> Adapted from an internal workflow guide prepared by **Mark Yang** and shared here with permission.  
> Paths, credentials, and institution-specific details have been removed or generalized.

## Scope

This document records the external workflow used for:

1. whole-slide image patching with CLAM,
2. feature extraction with a pretrained UNI2-h encoder,
3. HoVer-Net nuclei segmentation and classification, and
4. optional integration of CLAM attention scores with HoVer-Net cell-composition summaries.

The CLAM, UNI2-h, and HoVer-Net source code and pretrained weights are **not redistributed** in this repository.

---

# CLAM: whole-slide image preprocessing and feature extraction

## 1. Activate the CLAM environment

```bash
conda activate <clam_environment>
```

## 2. Prepare slide and output directories

```bash
mkdir -p <data_directory>
mkdir -p <results_directory>
```

Populate `<data_directory>` with symbolic links to the required `.svs` slides:

```bash
cd <data_directory>
ln -s <path_to_slides>/*.svs .
```

## 3. Segment tissue and extract patches

Run CLAM's `create_patches_fp.py`:

```bash
python <path_to_clam>/create_patches_fp.py   --source <data_directory>   --save_dir <results_directory>   --patch_size 256   --seg   --patch   --stitch
```

Generate a slide list for feature extraction:

```bash
ls <results_directory>/patches/ | sed 's/\.h5$//' > <results_directory>/h5_list.csv
```

## 4. Prepare the UNI2-h encoder

Access to UNI2-h requires registration with the model provider.

Store the model checkpoint path in an environment variable:

```bash
export UNI_CKPT_PATH=<path_to_uni2_h_checkpoint>
```

Do **not** hard-code or commit access tokens. Authenticate using an environment variable or the relevant command-line login flow.

The CLAM encoder builder must be configured to load the UNI2-h architecture and checkpoint. In the original institutional setup, the CLAM model name remained `uni_v1` even though UNI2-h weights were used.

## 5. Extract slide-level patch features

```bash
python <path_to_clam>/extract_features_fp.py   --data_h5_dir <results_directory>   --data_slide_dir <data_directory>   --csv_path <results_directory>/h5_list.csv   --feat_dir <feature_directory>   --batch_size <batch_size>   --slide_ext .svs   --model_name uni_v1
```

The resulting feature files are large and are not included in this repository.

---

# HoVer-Net: nuclei segmentation and classification

## 1. Activate the HoVer-Net environment

```bash
conda activate <hovernet_environment>
```

## 2. Obtain pretrained weights

The original analysis used a pretrained HoVer-Net checkpoint compatible with nuclei segmentation and classification.

Store model weights outside the repository:

```text
<path_to_hovernet_weights>
```

## 3. Define nuclei classes

The class mapping used in the analysis was:

```json
{
  "0": ["unlabeled",  [0, 0, 0]],
  "1": ["epithelial", [255, 0, 0]],
  "2": ["lymphocyte", [0, 255, 0]],
  "3": ["macrophage", [0, 0, 255]],
  "4": ["neutrophil", [255, 255, 0]]
}
```

Save this configuration as `type_info.json`.

## 4. Run whole-slide inference

```bash
python <path_to_hovernet>   --gpu=<gpu_ids>   --nr_types=5   --batch_size=<batch_size>   --model_mode=fast   --type_info_path=<path_to_type_info.json>   --model_path=<path_to_hovernet_checkpoint>   --nr_inference_workers=<n_workers>   --nr_post_proc_workers=<n_workers>   wsi   --input_dir=<input_slide_directory>   --output_dir=<output_directory>   --save_thumb   --save_mask   --proc_mag=<magnification>
```

The magnification and checkpoint must match the source slides and model configuration.

---

# Combined CLAM attention and HoVer-Net analysis

## 1. Extract attention scores

```bash
python <path_to_project>/get_attention.py
```

This step exports tile-level CLAM attention scores.

## 2. Summarize cell composition by attention level

```bash
python <path_to_project>/generate_cell_composition_summary.py
```

The original workflow grouped attention scores into low, medium, and high categories and calculated cell counts and percentages for each category.

## 3. Generate visual summaries

```bash
python <path_to_project>/plot_cell_composition_analysis.py
```

Potential outputs include:

- stacked bar charts of cell composition,
- box plots of cell counts,
- box plots of cell percentages, and
- scatter plots of attention score versus total cell count.

---

## Reproducibility note

This document records the institutional workflow used in the thesis. Exact reproduction depends on:

- the external CLAM and HoVer-Net codebases,
- the corresponding conda environments,
- pretrained model weights,
- local modifications to the CLAM encoder builder,
- institutional GPU infrastructure, and
- source whole-slide images that are not redistributed here.

Environment files and exact package versions can be added separately once recovered from the original compute environment.
