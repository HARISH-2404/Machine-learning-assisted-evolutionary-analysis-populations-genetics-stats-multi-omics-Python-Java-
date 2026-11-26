"""Shared validation, reproducibility, logging and output helpers."""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd

def get_logger(name: str) -> logging.Logger:
    """Return a consistently configured module logger."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger

def rng_from(seed: int | None = 20261127) -> np.random.Generator:
    """Create an isolated NumPy random generator."""
    return np.random.default_rng(seed)

def ensure_dir(path: str | Path) -> Path:
    """Create an output directory if necessary."""
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return target

def write_csv(frame: pd.DataFrame, path: str | Path) -> Path:
    """Write a DataFrame as CSV and return the resulting path."""
    target = Path(path)
    ensure_dir(target.parent)
    frame.to_csv(target, index=False)
    return target

def validate_diploid_dosage(frame: pd.DataFrame, id_column: str = "genotype_id") -> list[str]:
    """Validate a 0/1/2 dosage matrix; NaN is permitted as missing data."""
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Genotype input must be a pandas DataFrame.")
    if id_column not in frame.columns:
        raise ValueError(f"Missing required identifier column: {id_column}")
    if frame[id_column].isna().any() or frame[id_column].duplicated().any():
        raise ValueError(f"{id_column} must contain unique, non-missing IDs.")
    markers = [c for c in frame.columns if c != id_column]
    if not markers:
        raise ValueError("No marker columns found.")
    try:
        values = frame[markers].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Marker columns must contain numeric dosage values.") from exc
    observed = values[np.isfinite(values)]
    if not np.isin(observed, [0, 1, 2]).all():
        raise ValueError("Observed genotype dosage values must be 0, 1 or 2; missing values may be NaN.")
    return markers

def validate_probability(value: float, name: str) -> None:
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be between 0 and 1.")

def summarize_dosage(genotype: pd.DataFrame, id_column: str = "genotype_id") -> pd.DataFrame:
    """Calculate marker call rate, allele frequency, MAF and heterozygosity."""
    markers = validate_diploid_dosage(genotype, id_column)
    x = genotype[markers].astype(float)
    called = x.notna().sum(axis=0)
    af = x.sum(axis=0, skipna=True) / (2 * called.replace(0, np.nan))
    return pd.DataFrame({
        "marker_id": markers,
        "call_rate": called.to_numpy() / len(genotype),
        "alternate_allele_frequency": af.to_numpy(),
        "minor_allele_frequency": np.minimum(af, 1-af).to_numpy(),
        "observed_heterozygosity": (x == 1).sum(axis=0).to_numpy() / called.replace(0, np.nan).to_numpy()
    })
