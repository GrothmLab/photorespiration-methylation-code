#!/usr/bin/env bash
# Download and unpack Zenodo source-data archives into the expected R subfolders.
# Usage: edit RECID and FILENAMES below, then run:  bash download_source_data.sh
set -euo pipefail

if ! command -v curl >/dev/null 2>&1; then
  echo "Error: curl is required" >&2
  exit 1
fi
if ! command -v unzip >/dev/null 2>&1; then
  echo "Error: unzip is required" >&2
  exit 1
fi

RECID="18135710"  # TODO: set your Zenodo record ID
BASE_URL="https://zenodo.org/records/${RECID}/files"

# TODO: set your exact Zenodo file names
LD_SD_ZIP="LD_SD_source-data.zip"
HCO2_ZIP="hCO2_cCO2_source-data.zip"
MT_ZIP="m_t_mt_wt_source-data.zip"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

download_and_place() {
  local zipname="$1"
  local target_dir="$2"

  mkdir -p "$target_dir"
  echo "Downloading ${zipname} ..."
  curl -L -o "$TMPDIR/$zipname" "${BASE_URL}/${zipname}?download=1"

  echo "Unzipping ${zipname} ..."
  # Unzip into a temp extract dir to avoid wrapper-folder issues
  local extract_dir="$TMPDIR/extract_${zipname%.zip}"
  mkdir -p "$extract_dir"
  unzip -q "$TMPDIR/$zipname" -d "$extract_dir"

  # If the zip contains a single top-level folder, flatten it
  shopt -s nullglob dotglob
  local entries=("$extract_dir"/*)
  if [ ${#entries[@]} -eq 1 ] && [ -d "${entries[0]}" ]; then
    echo "Flattening single wrapper folder in ${zipname}"
    extract_dir="${entries[0]}"
  fi

  echo "Moving contents into ${target_dir} ..."
  shopt -s dotglob nullglob
  mv "$extract_dir"/* "$target_dir"/ || true
}

download_and_place "$LD_SD_ZIP" "R/WGBS_analysis_LD_SD/source-data"
download_and_place "$HCO2_ZIP" "R/WGBS_analysis_hCO2_cCO2/source-data"
download_and_place "$MT_ZIP" "R/WGBS_analysis_m_t_mt_wt/source-data"

echo "All done."
