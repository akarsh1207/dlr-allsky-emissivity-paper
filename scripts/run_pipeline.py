#!/usr/bin/env python3
"""Central runner for reader-facing all-sky DLR analysis workflows.

Usage example:
	python scripts/run_pipeline.py --h5 output.h5 --output outputs/plots
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


STEFAN_BOLTZMANN = 5.670374419e-8
STATIONS = ["BON", "DRA", "FPK", "GWN", "PSU", "SXF", "TBL"]


def _prepare_station(df: pd.DataFrame) -> pd.DataFrame:
	"""Apply core quality filters and compute derived fields used in plots/metrics."""
	df = df.copy()
	df.index = pd.to_datetime(df.index)
	df = df.replace(-9999.9, np.nan)

	req = ["sza", "ghi_m", "dhi_m", "dni_m", "ghi_c", "dni_c", "dlw", "temp", "rh"]
	df = df.dropna(subset=req)

	# Core physical and quality filters from manuscript workflow.
	mask = (
		(df["sza"] < 72.5)
		& (df["ghi_m"] > 0)
		& (df["dhi_m"] > 0)
		& (df["dni_m"] > 0)
		& ((df["ghi_m"] / df["ghi_c"]) >= 0.1)
		& ((df["ghi_m"] / df["ghi_c"]) <= 1.5)
		& ((df["dhi_m"] / df["ghi_m"]) <= 1)
		& (df["dlw"] > 0)
		& (df["temp"] >= -80)
		& (df["temp"] <= 90)
	)

	# If QC columns exist, enforce 0 for critical channels.
	for qc_col in ["qc_direct_n", "qc_dwsolar", "qc_diffuse", "qc_dwir", "qc_temp", "qc_rh", "qc_pressure"]:
		if qc_col in df.columns:
			mask &= df[qc_col] == 0

	df = df.loc[mask].copy()

	df["k_t"] = df["ghi_m"] / df["ghi_c"]
	df["k_d"] = df["dhi_m"] / df["ghi_m"]

	e_s = 6.112 * np.exp(17.625 * df["temp"] / (df["temp"] - 30.11 + 273.15))
	pw_hpa = e_s * df["rh"] / 100.0
	df = df.loc[pw_hpa > 0].copy()
	pw_hpa = pw_hpa.loc[df.index]

	e_clear = 0.6 + 1.652 * np.sqrt(pw_hpa / 1013.25)
	gamma = 0.585 * np.power(df["k_d"], 1.748)
	e_pred = (1 - gamma) * e_clear + gamma
	e_obs = df["dlw"] / (STEFAN_BOLTZMANN * (df["temp"] + 273.15) ** 4)

	out = pd.DataFrame(index=df.index)
	out["k_d"] = df["k_d"]
	out["gamma"] = gamma
	out["e_obs"] = e_obs
	out["e_pred"] = e_pred
	out = out.replace([np.inf, -np.inf], np.nan).dropna()
	out = out[(out["gamma"] <= 1) & (out["k_d"] <= 1)]
	return out


def _rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
	return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def _r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
	ss_res = np.sum((y_true - y_pred) ** 2)
	ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
	if ss_tot == 0:
		return float("nan")
	return float(1 - ss_res / ss_tot)


def _station_metrics(df: pd.DataFrame) -> dict[str, float]:
	y = df["e_obs"].to_numpy()
	yhat = df["e_pred"].to_numpy()
	rmse = _rmse(y, yhat)
	rrmse = float(np.sqrt(np.mean(((y - yhat) / y) ** 2)))
	mbe = float(np.mean(y - yhat))
	return {
		"count": int(len(df)),
		"mbe": mbe,
		"rmse": rmse,
		"rrmse": rrmse,
		"r2": _r2(y, yhat),
	}


def _plot_kd_gamma(df_all: pd.DataFrame, out_dir: Path) -> None:
	fig, ax = plt.subplots(figsize=(7, 5), dpi=180)
	hb = ax.hexbin(df_all["k_d"], df_all["gamma"], gridsize=80, cmap="viridis", mincnt=1)
	ax.set_xlabel("k_d")
	ax.set_ylabel("gamma")
	ax.set_title("Cloud Factor vs Diffuse Fraction")
	cb = fig.colorbar(hb, ax=ax)
	cb.set_label("counts")
	fig.tight_layout()
	fig.savefig(out_dir / "kd_vs_gamma_hexbin.png")
	plt.close(fig)


def _plot_obs_pred(df_all: pd.DataFrame, out_dir: Path) -> None:
	fig, ax = plt.subplots(figsize=(6, 6), dpi=180)
	ax.scatter(df_all["e_obs"], df_all["e_pred"], s=2, alpha=0.2)
	mn = min(df_all["e_obs"].min(), df_all["e_pred"].min())
	mx = max(df_all["e_obs"].max(), df_all["e_pred"].max())
	ax.plot([mn, mx], [mn, mx], "k--", lw=1)
	ax.set_xlabel("Observed emissivity")
	ax.set_ylabel("Predicted emissivity")
	ax.set_title("Observed vs Predicted Emissivity")
	fig.tight_layout()
	fig.savefig(out_dir / "observed_vs_predicted_emissivity.png")
	plt.close(fig)


def run(h5_path: Path, out_dir: Path, stations: Iterable[str]) -> None:
	out_dir.mkdir(parents=True, exist_ok=True)
	summaries: list[dict[str, float | str]] = []
	merged = []

	with pd.HDFStore(h5_path, mode="r") as store:
		available = set(k.strip("/") for k in store.keys())
		for st in stations:
			if st not in available:
				continue
			raw = store[st]
			prepared = _prepare_station(raw)
			if prepared.empty:
				continue
			m = _station_metrics(prepared)
			m["station"] = st
			summaries.append(m)
			prepared["station"] = st
			merged.append(prepared)

	if not merged:
		raise RuntimeError("No valid data remained after filtering. Check input file and filters.")

	df_all = pd.concat(merged, axis=0)
	_plot_kd_gamma(df_all, out_dir)
	_plot_obs_pred(df_all, out_dir)

	summary_df = pd.DataFrame(summaries).sort_values("station")
	summary_df.to_csv(out_dir / "emissivity_metrics_by_station.csv", index=False)

	overall = _station_metrics(df_all)
	overall["station"] = "ALL"
	pd.DataFrame([overall]).to_csv(out_dir / "emissivity_metrics_overall.csv", index=False)

	print(f"Saved outputs to: {out_dir}")
	print(f"Rows used (all stations): {len(df_all):,}")


def main() -> None:
	parser = argparse.ArgumentParser(description="Run reader-facing all-sky DLR analysis pipeline")
	parser.add_argument("--h5", required=True, type=Path, help="Path to SURFRAD HDF5 file (e.g., output.h5)")
	parser.add_argument(
		"--output",
		type=Path,
		default=Path("outputs/plots"),
		help="Output directory for plots/metrics",
	)
	parser.add_argument(
		"--stations",
		nargs="+",
		default=STATIONS,
		help="Stations to process (default: BON DRA FPK GWN PSU SXF TBL)",
	)
	args = parser.parse_args()
	run(args.h5, args.output, args.stations)


if __name__ == "__main__":
	main()
