# SURFRAD All-Sky DLR Dataset Package (v1.0)

This folder is a publication-ready package scaffold for Zenodo deposition.

## Recommended primary dataset filename

Rename `output.h5` to:

`surfrad_allsky_dlr_dataset_2010_2023_v1.0.h5`

Reasoning:
- Includes source network (`surfrad`)
- Includes purpose (`allsky_dlr_dataset`)
- Includes time coverage (`2010_2023` based on file index)
- Includes explicit version (`v1.0`)

## Files to include in Zenodo upload

Required:
- `surfrad_allsky_dlr_dataset_2010_2023_v1.0.h5` (renamed from `output.h5`)
- `docs/DATA_DICTIONARY.md`
- `docs/METHODS_SUMMARY.md`
- `metadata/zenodo_metadata_template.json`
- `scripts/quickstart_loader.py`
- `scripts/validate_hdf5_structure.py`

Optional:
- `all_precipitation_data.h5` (if required for reproducing rainy-day exclusion)
- `SHA256SUMS.txt` (generated with command shown below)

## Confirmed dataset structure from file inspection

- HDF5 root keys: BON, DRA, FPK, GWN, PSU, SXF, TBL
- Number of station tables: 7
- Columns per table: 51
- Datetime index with station-specific ranges from approximately 2010-01-01 to 2023-12-31
- Approximate file size: 14 GB

## Create checksum file

From this folder after placing the final HDF5 file:

```bash
shasum -a 256 surfrad_allsky_dlr_dataset_2010_2023_v1.0.h5 > SHA256SUMS.txt
```

## Final pre-upload checklist

1. Fill creator ORCIDs and add related identifiers (manuscript DOI/preprint and GitHub/software DOI) in `metadata/zenodo_metadata_template.json`
2. Confirm creator names, affiliations, ORCIDs, and funding
3. Set `access_right` to `embargoed` and set final embargo release date
4. Confirm DOI from Zenodo in manuscript data availability section
5. Upload files directly (do not rely only on a zip archive)
