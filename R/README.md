# R scripts for special BS-seq visualizations

This folder contains **R scripts for special visualizations** (e.g., clustering, heatmaps, genomic overlaps, metaplots, correlations) used in the associated study. 
Scripts assume the **working directory is the script's own folder** and use **relative paths** to input files.

## TL;DR

- **Code here; source data on Zenodo (DOI: 10.5281/zenodo.18135710)**.
- Generate or download and extract source data (e.g. using the hepler script).
- Run scripts in **RStudio** (or **R.app/RGui**) with the script’s folder as the working directory.

## Using Zenodo source data with this repository

**Source data:** Large computational source data needed to run these scripts were archived on **Zenodo** (DOI: 10.5281/zenodo.18135710).
The GitHub repository contains code only.
The R scripts expect source-data/ inside each analysis folder under this repo:
```
R/WGBS_analysis_m_t_mt_wt/source-data/...
R/WGBS_analysis_LD_SD/source-data/...
R/WGBS_analysis_hCO2_cCO2/source-data/...
```
**Download** and extract the Zenodo archive (DOI: 10.5281/zenodo.18135710) and place the source data files into the expected `source-data/` folders, or create **symlinks** from each analysis folder to the corresponding data (link the `source-data` folder, not individual files).

- Option A) Use the following **command line (Linux/macOS)** to copy the three source-data/ folders into place:

```bash
# from the repo root
curl -L -o source-data.tar.gz "https://zenodo.org/records/18135710/files/photorespiration_methylation_zenodo_source_data.tar.gz?download=1"
tar -xzf source-data.tar.gz

mkdir -p R/WGBS_analysis_LD_SD/source-data R/WGBS_analysis_hCO2_cCO2/source-data R/WGBS_analysis_m_t_mt_wt/source-data

rsync -a photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_LD_SD/source-data/      R/WGBS_analysis_LD_SD/source-data/
rsync -a photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_hCO2_cCO2/source-data/  R/WGBS_analysis_hCO2_cCO2/source-data/
rsync -a photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_m_t_mt_wt/source-data/  R/WGBS_analysis_m_t_mt_wt/source-data/
```
- Option B) Use the included helper script **download_source_data.sh** that **downloads, extracts, and creates symlinks**:

```bash
# from the repository root (macOS/Linux, Git Bash, or WSL)
bash download_source_data.sh
```

What it does:
- Downloads and extracts the tarball into `photorespiration_methylation_zenodo_source_data/`
- Creates symlinks:
  - `R/WGBS_analysis_LD_SD/source-data  →  photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_LD_SD/source-data`
  - `R/WGBS_analysis_hCO2_cCO2/source-data  →  photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_hCO2_cCO2/source-data`
  - `R/WGBS_analysis_m_t_mt_wt/source-data  →  photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_m_t_mt_wt/source-data`

**Windows notes:** Use **Git Bash** or **WSL**. Native PowerShell can run `bash download_source_data.sh`; creating symlinks may require Administrator privileges.

## How to run (RStudio)

These scripts are intended to be run **inside RStudio** (not via `Rscript`).

1. Open the script in RStudio.

2. Set the working directory, e.g. by running 
   ```r
   if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
   setwd(dirname(rstudioapi::getSourceEditorContext()$path))
   }
   ```

3. Make sure the source data files exist, e.g. with
   ```r
   dir.exists("../source-data")
   ```

4. Source or run the script (e.g., *Run All*). Relative paths in the script assume the working directory is the script’s folder.

## Running outside RStudio (R.app / RGui)

You can also run scripts in the base R GUI:
```r
# from the R Console window
source("path/to/your_script.R", chdir = TRUE)
```

## Inputs

Refer to inline comments of each R script for the expected inputs.

## Reproducibility

Per-script package environments are recorded in `session_info.txt` files located in each R script parent folder. 
These files capture the R version, platform, and loaded package versions used to generate the figures.

## License & citation

- License: MIT (see `LICENSE` in this folder)
- Please cite the associated manuscript **and** the Zenodo archive for this code (https://doi.org/10.5281/zenodo.18137716).

## Associated manuscript

**Hankover V. _et al._** *Photorespiration is linked to DNA methylation by formate as a one-carbon source.* (in press, 2026). DOI: _TBA_.



