#!/usr/bin/env bash
set -euo pipefail

RECID="18135710"  # Zenodo record ID (data)
TARBALL="photorespiration_methylation_zenodo_source_data.tar.gz"
URL="https://zenodo.org/records/${RECID}/files/${TARBALL}?download=1"

echo "[i] Downloading ${TARBALL} from Zenodo (record ${RECID})"
curl -L -o "${TARBALL}" "${URL}"

echo "[i] Extracting tarball at repo root..."
tar -xzf "${TARBALL}"

MODE="${1:-symlink}"   # usage: bash download_source_data.sh [symlink|copy]

if [[ "${MODE}" == "symlink" ]]; then
  echo "[i] Creating symlinks to extracted source-data (no copying)"
  rm -rf R/WGBS_analysis_LD_SD/source-data R/WGBS_analysis_hCO2_cCO2/source-data R/WGBS_analysis_m_t_mt_wt/source-data
  ln -snf ../../photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_LD_SD/source-data      R/WGBS_analysis_LD_SD/source-data
  ln -snf ../../photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_hCO2_cCO2/source-data  R/WGBS_analysis_hCO2_cCO2/source-data
  ln -snf ../../photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_m_t_mt_wt/source-data  R/WGBS_analysis_m_t_mt_wt/source-data
else
  echo "[i] Copying source-data folders into place (may take time/space)"
  mkdir -p R/WGBS_analysis_LD_SD/source-data R/WGBS_analysis_hCO2_cCO2/source-data R/WGBS_analysis_m_t_mt_wt/source-data
  rsync -a photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_LD_SD/source-data/      R/WGBS_analysis_LD_SD/source-data/
  rsync -a photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_hCO2_cCO2/source-data/  R/WGBS_analysis_hCO2_cCO2/source-data/
  rsync -a photorespiration_methylation_zenodo_source_data/R/WGBS_analysis_m_t_mt_wt/source-data/  R/WGBS_analysis_m_t_mt_wt/source-data/
fi

echo "[✓] Done."