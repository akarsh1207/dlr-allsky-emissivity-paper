#!/usr/bin/env bash
set -euo pipefail

# Build a clean upload folder with renamed data file and documentation.
# Usage:
#   bash prepare_zenodo_bundle.sh /absolute/path/to/output.h5 [optional/path/to/all_precipitation_data.h5]

if [[ $# -lt 1 ]]; then
  echo "Usage: bash prepare_zenodo_bundle.sh /path/to/output.h5 [path/to/all_precipitation_data.h5]"
  exit 1
fi

OUTPUT_H5="$1"
PRECIP_H5="${2:-}"

if [[ ! -f "$OUTPUT_H5" ]]; then
  echo "Error: output file not found: $OUTPUT_H5"
  exit 1
fi

TARGET_DIR="zenodo_upload_ready_v1"
mkdir -p "$TARGET_DIR/docs" "$TARGET_DIR/metadata" "$TARGET_DIR/scripts"

cp "README_data.md" "$TARGET_DIR/"
cp "docs/DATA_DICTIONARY.md" "$TARGET_DIR/docs/"
cp "docs/METHODS_SUMMARY.md" "$TARGET_DIR/docs/"
cp "metadata/zenodo_metadata_template.json" "$TARGET_DIR/metadata/"
cp "metadata/zenodo_form_copy_paste.md" "$TARGET_DIR/metadata/"
cp "scripts/quickstart_loader.py" "$TARGET_DIR/scripts/"
cp "scripts/validate_hdf5_structure.py" "$TARGET_DIR/scripts/"

cp "$OUTPUT_H5" "$TARGET_DIR/surfrad_allsky_dlr_dataset_2010_2023_v1.0.h5"

if [[ -n "$PRECIP_H5" ]]; then
  if [[ -f "$PRECIP_H5" ]]; then
    cp "$PRECIP_H5" "$TARGET_DIR/all_precipitation_data.h5"
  else
    echo "Warning: precipitation file not found: $PRECIP_H5"
  fi
fi

(
  cd "$TARGET_DIR"
  shasum -a 256 surfrad_allsky_dlr_dataset_2010_2023_v1.0.h5 > SHA256SUMS.txt
)

echo "Created $TARGET_DIR with renamed dataset and support files."
