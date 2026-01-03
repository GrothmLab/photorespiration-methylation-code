# Code for “Photorespiration is linked to DNA methylation by formate as a one-carbon source” (in press)

This repository contains:
- `snakemake/` — BS-seq processing pipeline snapshot.
- `R/` — R scripts for analysis and visualizations of WGBS data.  
  See `R/README.md` for run instructions and data layout.

## Data & reproducibility
- **Source data** to run the R scripts are archived on **Zenodo**: https://doi.org/10.5281/zenodo.18135710  
  The GitHub repo contains code only.
- Each R script folder includes a `session_info.txt` with R and package versions.

## Quick start
- **R scripts**: open in RStudio (or R.app/RGui), set the working directory to the script’s folder, and run.  
  Scripts expect `source-data/` subfolders next to each analysis folder; see `R/README.md`.
- **Download data (CLI)**: from repo root, `bash download_source_data.sh` after filling your Zenodo record ID and filenames.

## Citation
Please cite:
1. **Associated manuscript**:  
   Hankover V. _et al._ *Photorespiration is linked to DNA methylation by formate as a one-carbon source.* (in press, 2026). DOI: _TBA_.
2. **This code/data archive** (Zenodo DOI): 10.5281/zenodo.18137716.  
   See also `CITATION.cff`.

## License
MIT (see license files). Third-party tools/scripts retain their original licenses.