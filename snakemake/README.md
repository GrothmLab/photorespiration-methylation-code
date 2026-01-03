# BS-seq Snakemake workflow

Snakemake workflow for bisulfite sequencing (BS-seq) processing and summary.

## Contents
- `workflow/Snakefile` — main workflow
- `config/config.yaml` — configuration
- `config/samples.template.tsv` — sample sheet template
- `tools` - helper scripts
- `README.md` — this file
- `versions.txt` - software version / requirements 

## What it does
-  Read trimming (Trimmomatic)
-  Alignment & deduplication (Bismark + Bowtie2)
-  Methylation extraction / reports (Bismark)
-  Methylkit input files
-  QC reports (FastQC, MultiQC)
-  Optional track conversion (e.g., bedGraph/bigWig; Bismark_to_wig.pl if used)

## Quickstart
1) Clone or download this repository.
2) Adjust paths in `config/config.yaml` and the sample sheet (`config/samples.template.tsv`).
3) Ensure conda or container support is available for Snakemake.

### Run
Example (local, with conda):
```bash
snakemake -s workflow/Snakefile --use-conda --cores 4 --configfile config/config.yaml
```

## Inputs
- Raw data paths defined in `config/config.yaml`
- Sample metadata in `config/samples.template.tsv` (headers included; replace dummy entries)
- Raw WGBS data from the associated study *Photorespiration is linked to DNA methylation by formate as a one-carbon source.* (in press, 2026) is deposited in the
  Gene Expression Omnibus (GEO) under the accession code GSE292915 https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE292915&format=file

## Outputs
- Workflow creates intermediate/output files under your configured work/results directories (see rules).

## Reproducibility
- Software versions are recorded ().

## Dependencies
This repository does **not** bundle third‑party tools. Install the software referenced in the Snakefile/rules yourself and ensure each executable is on your **PATH** (e.g., Trimmomatic, Bowtie2, Bismark, Samtools, FastQC, MultiQC, Snakemake). See `versions_final.csv` / `versions.txt` for versions used.

**Included helper scripts (provided for convenience; not auto‑referenced):**
- `Bismark_to_wig.pl`
- `filter_non_conversion_H` (CHH‑only variant of Bismark’s non‑conversion filter)

Copy these scripts to a directory on your **PATH** (e.g., `~/bin`) **or** update your rules to reference their full paths. They are **not** called from a local `tools/` folder by default.

## Notes & limitations
-  This snapshot is not a turnkey pipeline; environments are not pinned.
-  Replace any remaining personal paths (search for PLACEHOLDER_ or absolute paths).
-  If you prefer conda, add envs and run with --use-conda.

## License
Code in this repository is released under the MIT License (see `LICENSE`). Third‑party tools/scripts retain their original licenses.

## Citation
Please cite
1) Code snapshot (this GitHub release): https://doi.org/10.5281/zenodo.18137716
2) Associated manuscript: Hankover V. _et al._ *Photorespiration is linked to DNA methylation by formate as a one-carbon source.* (in press, 2026). DOI: _TBA_.

