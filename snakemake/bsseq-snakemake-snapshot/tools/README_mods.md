# Bismark `filter_non_conversion` modification — CHH-only filter

**Upstream tool:** Bismark (GPL-3.0-or-later)  
**Original file:** `filter_non_conversion`  
**Modified file:** `filter_non_conversion_H`  

## What changed
The original script flags/removes reads based on **non‑CG methylation** (i.e. **CHG _and_ CHH**) either by an absolute **`--threshold`** of methylated non‑CG calls or by a **`--percentage_cutoff`** (with `--minimum_count`).  
This modified version **restricts the counting to CHH context only**. In practice:
- **Original logic:** counts `H`/**`X`** as methylated non‑CG and `h`/**`x`** as unmethylated non‑CG (i.e. **CHH + CHG**).  
- **Modified logic:** counts only **`H`/`h`** (**CHH**) and **ignores `X`/`x`** (**CHG**).  
- Consequence: **Consecutive CHG methylation no longer triggers filtering**; only **consecutive CHH methylation** is considered for non‑conversion filtering (both in single‑end and paired‑end branches).

> Note: The help/report strings still mention “non‑CG (CHH or CHG)” from upstream. The functional change is the counter itself, which now targets **CHH‑only**.

## Why
To specifically remove likely non‑conversion artefacts arising from **CHH** context without penalizing legitimate or condition‑specific **CHG** methylation.

## Impact on results
- Reads with **runs of CHH methylation** are filtered as before.  
- Reads with **runs of CHG methylation** that would have been filtered by the original script are **retained** in this variant.  
- CpG methylation handling is unchanged. Output filenames/headers remain the upstream defaults.

## How to apply to upstream Bismark
1. Obtain the upstream Bismark release matching your environment (e.g., `v0.23.1`).  
2. Back up the original `filter_non_conversion`.  
3. Replace it with `filter_non_conversion_H`, **or** apply the accompanying unified diff.

A unified diff (`filter_non_conversion_H.patch`) from the original to the CHH‑only version is provided in this archive.

## Reproducibility notes
- Used in the BS‑seq analyses for the manuscript.  
- Both the **modified script** and the **patch** are included for transparency.  
- Toolchain example from logs: **Bismark v0.23.1** (see Zenodo entry for full versions).

## License
This modification and patch are distributed under **GPL‑3.0‑or‑later**, consistent with Bismark’s license. Upstream copyright remains with the original authors.
