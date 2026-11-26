"""Missingness filters and simple imputation helpers."""
import pandas as pd

def filter_missing_features(df, max_missing=0.2):
    if not 0 <= max_missing <= 1: raise ValueError("max_missing must be in [0,1]")
    keep=df.isna().mean(axis=0) <= max_missing
    return df.loc[:, keep].copy(), list(df.columns[~keep])

def impute_block(df, strategy="median"):
    """Simple per-feature imputation; fit this only on training folds in prediction tasks."""
    if strategy == "median": values=df.median(numeric_only=True)
    elif strategy == "mean": values=df.mean(numeric_only=True)
    elif strategy == "zero": values=pd.Series(0, index=df.columns)
    else: raise ValueError("strategy must be median, mean, or zero")
    return df.fillna(values), values
