# Methods Summary (for Zenodo record)

## Dataset purpose
This dataset supports estimation of all-sky atmospheric emissivity and downwelling longwave radiation (DLR) from shortwave and meteorological observations.

## Data source
- NOAA SURFRAD stations: BON, DRA, FPK, GWN, PSU, SXF, TBL
- Variables include radiometric observations, meteorology, quality-control flags, and clear-sky irradiance estimates.

## Temporal coverage
- Per file index inspection (exact by station):
	- BON: 2010-01-01 07:19:36 to 2023-12-31 16:43:03 (3,671,601 rows)
	- DRA: 2009-12-31 16:13:01 to 2023-12-31 16:12:28 (3,603,540 rows)
	- FPK: 2010-01-01 07:47:41 to 2023-12-31 16:14:08 (3,693,809 rows)
	- GWN: 2010-01-01 07:03:36 to 2023-12-31 16:59:03 (3,611,888 rows)
	- PSU: 2010-01-01 07:21:22 to 2023-12-31 16:40:49 (3,668,769 rows)
	- SXF: 2010-01-01 07:30:36 to 2023-12-31 16:31:03 (3,567,519 rows)
	- TBL: 2010-01-01 07:19:08 to 2023-12-31 16:42:35 (3,673,243 rows)
- Total rows in this release: 25,490,369.
- Note: this release covers data through 2023; if manuscript text mentions 2024 validation, that validation subset is not included in this specific file.

## Core processing (high level)
1. Quality control filtering using station-provided QC flags and physical plausibility thresholds.
2. Removal of rainy periods/dates based on precipitation criteria.
3. Computation of derived quantities used in modeling (for manuscript workflows).
4. Storage into station-wise tables in HDF5 format.

## Core processing (implemented thresholds)
- Solar-geometry/radiation filters:
	- sza < 72.5 degrees
	- ghi_m > 0, dhi_m > 0
	- 0.1 < ghi_m / ghi_c < 1.5
	- ghi_m < 1.2 * G_on * cos(sza)^1.2 + 50
	- dni_m < 0.95 * G_on * cos(sza)^0.2 + 10
	- k_d = dhi_m / ghi_m <= 1
	- local time >= 08:00
- Meteorology and QC filters:
	- dlw > 0
	- -80 <= temp <= 90 (degC)
	- e_s >= 0 and pw_hpa >= 0
	- required QC flags equal to 0: qc_direct_n, qc_dwsolar, qc_diffuse, qc_dwir, qc_temp, qc_rh, qc_pressure
- Derived-quantity filter:
	- r = (e_sky - e_clear_sky) / (1 - e_clear_sky) <= 1
- Cleanup:
	- remove inf/-inf and missing values after derived calculations

## HDF5 layout
- One top-level group per station.
- Main table per station under `<STATION>/table` with datetime index and 51 columns.

## Reproducibility notes
- Use scripts in `scripts/` for quick loading and structural validation.
- Full analysis/training scripts should be maintained in the companion public GitHub repository.

## Software environment used for this package
- Python: 3.12 (project virtual environment)
- Core libraries used in support scripts: pandas, PyTables/HDF5

## Final items to complete before publish
- Add final public GitHub repository URL in Zenodo related identifiers and manuscript Data Availability text.
- Add ORCID IDs for all creators in metadata.
