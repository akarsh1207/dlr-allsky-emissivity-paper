#!/usr/bin/env python3
"""Quick loader for station-wise tables in the SURFRAD HDF5 dataset."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Load and preview a station table from HDF5")
    parser.add_argument("h5_path", type=Path, help="Path to HDF5 dataset")
    parser.add_argument("--station", default="BON", help="Station key: BON, DRA, FPK, GWN, PSU, SXF, TBL")
    parser.add_argument("--rows", type=int, default=5, help="Number of rows to preview")
    args = parser.parse_args()

    station_key = f"/{args.station.upper()}"

    with pd.HDFStore(args.h5_path, mode="r") as store:
        keys = store.keys()
        if station_key not in keys:
            raise ValueError(f"Station {station_key} not found. Available keys: {keys}")
        df = store.select(station_key, start=0, stop=args.rows)

    print(f"Loaded station: {station_key}")
    print(f"Shape preview: {df.shape}")
    print(df.head(args.rows))


if __name__ == "__main__":
    main()
