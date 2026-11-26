"""Shared I/O, validation and reproducible preprocessing helpers."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

def save_table(frame: pd.DataFrame, path: str | Path) -> Path:
    target=Path(path)
    target.parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(target,index=False)
    return target

def require_columns(frame: pd.DataFrame, required, label="table"):
    missing=sorted(set(required)-set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")

def require_unique(frame: pd.DataFrame, column: str):
    require_columns(frame,[column])
    if frame[column].isna().any() or frame[column].duplicated().any():
        raise ValueError(f"{column} must contain unique non-null identifiers.")

def numeric_matrix(frame: pd.DataFrame, columns):
    try:
        matrix=frame[list(columns)].apply(pd.to_numeric,errors="raise").to_numpy(dtype=float)
    except (TypeError,ValueError) as exc:
        raise ValueError("Feature columns must contain numeric values.") from exc
    if np.isinf(matrix).any():
        raise ValueError("Infinite feature values are not allowed.")
    return matrix
