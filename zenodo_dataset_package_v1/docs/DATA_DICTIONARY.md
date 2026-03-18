# Data Dictionary

This dictionary corresponds to station tables stored in the HDF5 file `surfrad_allsky_dlr_dataset_2010_2023_v1.0.h5`.

## Station groups
- BON
- DRA
- FPK
- GWN
- PSU
- SXF
- TBL

## Table schema
Each station table contains 51 float columns and a datetime-like index.

| Column | Type | Full Range (all stations) | Unit | Description |
|---|---|---:|---|---|
| year | float64 | 2010 to 2023 | year | Calendar year |
| month | float64 | 1 to 12 | month | Calendar month |
| day | float64 | 1 to 31 | day | Day of month |
| jday | float64 | 1 to 366 | day | Julian day |
| hour | float64 | 0 to 23 | hour | Hour of day |
| min | float64 | 0 to 59 | minute | Minute of hour |
| dt | float64 | 0 to 23.983 | hour | Decimal local time of day (approximately hour + min/60, rounded) |
| sza | float64 | 10.81 to 91.0 | degree | Solar zenith angle |
| dni_m | float64 | -10.0 to 1196.8 | W m^-2 | Measured direct normal irradiance |
| ghi_m | float64 | -9999.9 to 1541.1 | W m^-2 | Measured global horizontal irradiance |
| uw_solar | float64 | -9999.9 to 1013.8 | W m^-2 | Upwelling shortwave irradiance |
| dhi_m | float64 | -19.0 to 1202.5 | W m^-2 | Measured diffuse horizontal irradiance |
| dlw | float64 | 118.2 to 546.6 | W m^-2 | Downwelling longwave irradiance |
| dw_casetemp | float64 | 233.86 to 321.2 | K | Downwelling pyrgeometer case temperature |
| dw_dometemp | float64 | 233.65 to 320.9 | K | Downwelling pyrgeometer dome temperature |
| uw_ir | float64 | -9999.9 to 917.4 | W m^-2 | Upwelling longwave irradiance |
| uw_casetemp | float64 | -9999.9 to 321.35 | K | Upwelling pyrgeometer case temperature |
| uw_dometemp | float64 | -9999.9 to 406.3 | K | Upwelling pyrgeometer dome temperature |
| uvb | float64 | -9999.9 to 354.2 | SURFRAD native | UVB channel as reported by SURFRAD source stream |
| par | float64 | -9999.9 to 645.9 | SURFRAD native | Photosynthetically active radiation channel as reported by SURFRAD source stream |
| netsolar | float64 | -9999.9 to 1361.1 | W m^-2 | Net shortwave irradiance |
| netir | float64 | -9999.9 to 0.0 | W m^-2 | Net longwave irradiance |
| totalnet | float64 | 0.0 to 0.0 | W m^-2 | Total net radiation |
| temp | float64 | -9999.9 to 174.4 | degC | Air temperature |
| rh | float64 | -9999.9 to 260.5 | % | Relative humidity |
| windspd | float64 | -9999.9 to 77.6 | m s^-1 | Wind speed |
| winddir | float64 | -9999.9 to 360.0 | degree | Wind direction |
| pressure | float64 | -9999.9 to 1031.9 | hPa | Surface pressure |
| qc_direct_n | float64 | 0.0 to 0.0 | flag | QC flag |
| qc_netsolar | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_netir | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_dwsolar | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_uwsolar | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_diffuse | float64 | 0.0 to 0.0 | flag | QC flag |
| qc_dwir | float64 | 0.0 to 0.0 | flag | QC flag |
| qc_dwcasetemp | float64 | 0.0 to 0.0 | flag | QC flag |
| qc_dwdometemp | float64 | 0.0 to 0.0 | flag | QC flag |
| qc_uwir | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_uwcasetemp | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_uwdometemp | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_uvb | float64 | 0.0 to 2.0 | flag | QC flag |
| qc_par | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_totalnet | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_temp | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_rh | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_windspd | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_winddir | float64 | 0.0 to 1.0 | flag | QC flag |
| qc_pressure | float64 | 0.0 to 0.0 | flag | QC flag |
| ghi_c | float64 | 0.0 to 1051.526168946726 | W m^-2 | Clear-sky GHI estimate |
| dhi_c | float64 | 0.0 to 192.24128799796756 | W m^-2 | Clear-sky DHI estimate |
| dni_c | float64 | 0.0 to 980.4605347103717 | W m^-2 | Clear-sky DNI estimate |

## Missing/fill values
- The value `-9999.9` appears in multiple variables and should be treated as missing/fill where appropriate.
- In downstream workflows, nonzero QC flags and physically impossible values are filtered.

## QC semantics used in manuscript workflows
- QC filters are applied with equality to zero for critical variables (for example qc_direct_n, qc_dwsolar, qc_diffuse, qc_dwir, qc_temp, qc_rh, qc_pressure).
- Operational interpretation used in scripts: 0 = accepted quality, nonzero = flagged.

## Notes
- Ranges above are raw data ranges across all stations and include fill values and occasional outliers.
- For modeling, apply the documented physical and QC filters before computing derived variables.
