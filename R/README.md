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

**Download** and extract the Zenodo archive (DOI: 10.5281/zenodo.18135710) and place the source data files into the expected `source-data/` folders, or create **symlinks** from each analysis folder to the corresponding data (link the `source-data` folder, not individual files). If you prefer command-line, you can fetch the Zenodo zips and extract them into the expected source-data/ folders:

- Option A) Use the following **command line (Linux/macOS)** — replace the URLs with the actual Zenodo file links:

```bash
# from the repo root
curl -L -o LD_SD.zip     "https://zenodo.org/records/18135710/files/LD_SD_source-data.zip?download=1" \
  && unzip -q LD_SD.zip -d R/WGBS_analysis_LD_SD/source-data \
  && curl -L -o hCO2_cCO2.zip "https://zenodo.org/records/18135710/files/hCO2_cCO2_source-data.zip?download=1" \
  && unzip -q hCO2_cCO2.zip -d R/WGBS_analysis_hCO2_cCO2/source-data \
  && curl -L -o m_t_mt_wt.zip "https://zenodo.org/records/18135710/files/m_t_mt_wt_source-data.zip?download=1" \
  && unzip -q m_t_mt_wt.zip -d R/WGBS_analysis_m_t_mt_wt/source-data
```
- Option B) using the included **Helper script:** run `bash download_source_data.sh` from the repo root after filling the `RECID` and file names.

Source data can be reproduced from raw data, available on Gene Expression Omnibus (GEO) under the accession codes GSE292915 and GSE292917.

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



