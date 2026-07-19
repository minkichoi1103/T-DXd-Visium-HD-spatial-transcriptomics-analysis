# T-DXd Visium HD spatial-transcriptomics analysis

Analysis pipeline for the study
***"Spatial transcriptomics uncovers pharmacokinetic barriers and tumor-intrinsic determinants of resistance to trastuzumab deruxtecan in breast cancer."***

The repository is organized as an ordered, step-by-step workflow. Notebooks and scripts
are numbered `1 … 11` and grouped into numbered folders by analysis stage. Running them in
order reproduces the intermediate tables and the final statistical analysis.

> **Data is not included.** This repository contains **code only**. The raw Visium HD
> `.h5ad` files, the clinical metadata file, and all generated intermediates are expected
> under `data/` (see [Data layout](#data-layout)); figures are written to `figures/`.
> The **raw Visium HD data** are available from the Zenodo community:
> <https://zenodo.org/communities/tdxd-visiumhd/records>. The clinical metadata is available
> as **Supplementary Data 1** of the paper — save it as `data/clinical.csv` (used by steps 1 and 4).

## Pipeline

Every step reads its inputs from `../data/` and writes outputs to `../data/` (tables /
`.h5ad`) or `../figures/` (plots), relative to the step's own folder. **Run each notebook
or script from the folder it lives in.**

| # | Folder / file | What it does | Reads → Writes (under `data/`) |
|---|---|---|---|
| 1 | `1.Quality control and filtering/`<br>`1_Quality control and filtering.ipynb` | Load per-sample Visium HD, split by core, compute the fraction of bins passing a range of total-count cutoffs, and plot QC curves. | `h5ad/`, `clinical.csv` → `figures/` |
| 2 | `2.Tissue fragment and segmentation/`<br>`2_Tissue fragment and segmentation.ipynb` | Split each median-filtered slide into tissue fragments (per slide / core / subcluster) and save one `.h5ad` per fragment. | `h5ad_segmented_mf/` → `h5ad_batch/` |
| 3 | `3.Segmentation of tumor cells and vessel area …/`<br>`3_median_filtering.ipynb`, `processing.py` | NaN-aware median filtering of cell-type abundance maps. `processing.py` provides `find_cell_areas` (Otsu thresholding → cell / edge / background areas) and `select_specific_cell`. | `h5ad_segmented/` → `h5ad_segmented_mf/` |
| 4 | `4.Spatial autocorrelation of HER2 …/`<br>`4_Cancer_area_detection_and_Spatial_autocorrelation.ipynb` | Detect cancer areas, build cancer pseudobulk, compute Geary's C and Moran's I of `ERBB2` per fragment, and merge clinical metadata. | `h5ad_batch/`, `clinical.csv` → `pseudobulk_cancer.csv`, `results_df_spatial.csv`, `h5ad_cancer/` |
| 5 | `4.Spatial autocorrelation of HER2 …/`<br>`5_Calculate the cancer vessel diatance.ipynb` | Compute cancer-to-vessel distances (KD-tree) per fragment and merge them into the table. | `results_df_spatial.csv` → `results_df_spatial_vessel.csv`, `figures/` |
| 6 | `5.Detection of HER2 and cathepsin …/`<br>`6_Tissue fragment and segmentation.ipynb` | Fragment classification for the 55 µm pseudo-Visium grid used by the PBPK model. | `h5ad_pbpk_55um/` → `h5ad_batch_55um/` |
| 7 | `5.Detection of HER2 and cathepsin …/`<br>`7_stopover_pseudovisium.ipynb` | STopover colocalization of `ERBB2` with cathepsins (CTSB / CTSL, the ADC linker enzymes) per fragment; merge into the table. | `results_df_spatial_vessel.csv` → `results_df_spatial_vessel_linkerenz.csv` |
| 8 | `6.Mathematical modeling of ADC payload …/`<br>`Generate_55um_pseudovisium.py`, `8_PBPK_55um.ipynb` | Generate a 55 µm grid pseudo-Visium representation, then run a PBPK ADC simulation of payload distribution (`pbpk_st.PBPKSimulation_ADC`). | `h5ad_segmented_mf/` → `h5ad_55um_merged/` → `h5ad_pbpk_55um_subcluster_endo/` |
| 9 | `7.Differential expression … and Pathway activity …/`<br>`9_pseudobulk_cancer_decoupler.R` | limma differential expression (resistant vs. sensitive) on cancer pseudobulk and decoupleR pathway-activity inference; bubble plots. | `pseudobulk_cancer.csv`, `results_df_spatial_vessel_linkerenz.csv` → `pathway_activity_decoupler.csv`, `DEG_results_*.csv`, `Supplementary_table2.csv`, `figures/` |
| 10 | `8.Ligand-receptor analysis …/`<br>`10_commot.ipynb`, `10_commat_rra_result.R` | COMMOT spatial ligand–receptor communication in peritumoral areas; robust rank aggregation (RRA) and plotting in R. | `h5ad_fragment/`, `results_df_spatial_vessel_linkerenz.csv` → `commot_ranked_result.csv`, `figures/` |
| 11 | `9.Statistical analysis/`<br>`11_statistical_anslysis.ipynb` | GEE (exchangeable working correlation) of cancer–vessel distance vs. clinical response (sensitive / resistant), stratified by HER2 status; violin plots. | `results_df_spatial_vessel_linkerenz.csv` → `figures/` |

### Metadata table flow

The per-fragment metadata table is built up across steps and consumed by steps 9, 10, 11:

```
4  ──► results_df_spatial.csv
5  ──► results_df_spatial_vessel.csv
7  ──► results_df_spatial_vessel_linkerenz.csv  ──►  9 · 10 · 11
```

Cancer pseudobulk from step 4 (`pseudobulk_cancer.csv`) feeds the step-9 DEG analysis.
Ranked LR results from step 10's notebook (`commot_ranked_result.csv`) feed its R script.

## Data layout

Create a `data/` folder next to the step folders and populate it with your own inputs.
Both the user-supplied raw inputs and the pipeline-generated intermediates live here:

```
T_Dxd_VisiumHD/
├── data/                 # inputs + generated intermediates (git-ignored)
│   ├── h5ad/                     # per-sample Visium HD (input, step 1)
│   ├── h5ad_segmented/           # segmented (input, step 3)
│   ├── h5ad_segmented_mf/        # median-filtered (steps 3→2/8)
│   ├── h5ad_batch/               # per-fragment (steps 2→4/5/7)
│   ├── h5ad_batch_55um/          # 55 µm fragments (steps 6→7)
│   ├── h5ad_pbpk_55um/           # PBPK grid input (input, step 6)
│   ├── h5ad_cancer/ , h5ad_fragment/
│   ├── clinical.csv              # clinical metadata — paper Supplementary Data 1 (inputs, steps 1 & 4)
│   └── results_df_*.csv , pseudobulk_cancer.csv , ...  # generated
├── figures/              # all plot outputs (git-ignored)
└── <step folders>/
```

## Requirements

**Python** (see `requirements.txt`): scanpy, squidpy, anndata, pandas, numpy, scipy,
scikit-image, matplotlib, seaborn, statsmodels, tacco, commot, and
[STopover](https://github.com/bsungwoo/STopover) (image-based ST topological colocalization).

**R** (steps 9 and 10): `limma`, `decoupleR`, the `tidyverse` (dplyr, readr, tibble,
forcats, ggplot2), `ggsci`, and `RobustRankAggreg` (for the RRA step).

### Proprietary tools (Portrai Inc.)

Two steps rely on internal Portrai Inc. software that is not distributed here; both are
described in the paper's Methods:

- **`pbpk_st` (`PBPKSimulation_ADC`)** — part of Portrai Inc.'s **TME-PK** program, used for the
  ADC payload-distribution simulation (step 8). See Methods,
  *"Mathematical modeling of antibody-drug conjugates payload distribution."*
- **`curvsplit`** — an internal Portrai Inc. algorithm that selects and separates the aligned
  regions of core-needle-biopsy specimens to define tissue fragments (step 2). See Methods,
  *"Tissue fragment and segmentation."*

## Notes / caveats

- **Run each file from its own folder** so the `../data/` and `../figures/` relative
  paths resolve correctly.
- Notebooks 4, 5, and 7 import the shared `processing.py` module (kept in the
  `3.Segmentation of tumor cells and vessel area …` folder) via a `sys.path` line at the
  top of each notebook.
- Notebook cell outputs have been cleared; re-run to regenerate figures and tables.

## License

License **to be decided** Until a license is added, this code is **all rights reserved** by the authors.
