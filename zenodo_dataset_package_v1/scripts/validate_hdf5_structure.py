#!/usr/bin/env python3
"""Validate expected schema for SURFRAD HDF5 release file."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

EXPECTED_STATIONS = ["/BON", "/DRA", "/FPK", "/GWN", "/PSU", "/SXF", "/TBL"]
EXPECTED_COLS = 51


def _index_is_datetime_like(index: pd.Index) -> bool:
    """Return True if index is DatetimeIndex or parseable timestamp values."""
    if isinstance(index, pd.DatetimeIndex):
        return True
    try:
        parsed = pd.to_datetime(index, errors="coerce")
        return parsed.notna().all()
    except Exception:
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate HDF5 file structure")
    parser.add_argument("h5_path", type=Path, help="Path to HDF5 dataset")
    args = parser.parse_args()

    errors: list[str] = []

    with pd.HDFStore(args.h5_path, mode="r") as store:
        keys = sorted(store.keys())
        if keys != sorted(EXPECTED_STATIONS):
            errors.append(f"Station keys mismatch. Found {keys}, expected {EXPECTED_STATIONS}")

        for key in EXPECTED_STATIONS:
            if key not in keys:
                continue
            storer = store.get_storer(key)
            nrows = storer.nrows
            if nrows <= 0:
                errors.append(f"{key}: empty table")
                continue
            sample = store.select(key, start=0, stop=1)
            if sample.shape[1] != EXPECTED_COLS:
                errors.append(f"{key}: expected {EXPECTED_COLS} columns, found {sample.shape[1]}")
            if not _index_is_datetime_like(sample.index):
                errors.append(f"{key}: index is not datetime-like")

    if errors:
        print("VALIDATION FAILED")
        for err in errors:
            print("-", err)
        raise SystemExit(1)

    print("VALIDATION PASSED")


if __name__ == "__main__":
    main()
