# All-Sky DLR Estimation Scripts

Reader-facing code repository for the paper:
"Estimation of Longwave All-sky Emissivity using Shortwave and Meteorological Measurements".

## Repository Layout

- `scripts/run_pipeline.py`: central runner for reproducible filtering, metrics, and plotting
- `scripts/legacy/`: original research scripts used during model development
- `src/dlr/`: reusable DLR estimation module
- `notebooks/`: analysis notebooks
- `outputs/plots/`: default destination for generated plots and metric CSV files

## Input Data

Place or reference your SURFRAD HDF5 file (for example `output.h5`) locally.
Large data files are intentionally not versioned in this repository.

## Quick Start

1. Create and activate a Python environment.
2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Run the central pipeline:

```bash
python scripts/run_pipeline.py --h5 /path/to/output.h5 --output outputs/plots
```

## Pipeline Outputs

The central runner writes the following by default:

- `outputs/plots/kd_vs_gamma_hexbin.png`
- `outputs/plots/observed_vs_predicted_emissivity.png`
- `outputs/plots/emissivity_metrics_by_station.csv`
- `outputs/plots/emissivity_metrics_overall.csv`

## Citation

If you use this code, please cite the associated paper and the Zenodo dataset record.
