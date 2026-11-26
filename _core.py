"""Internal shared validation and calculations."""
from __future__ import annotations
import numpy as np
import pandas as pd

def validate_dosage(df, id_col="genotype_id", population_col=None):
    if not isinstance(df,pd.DataFrame): raise TypeError("Expected pandas DataFrame.")
    for c in [id_col]+([population_col] if population_col else []):
        if c not in df.columns: raise ValueError(f"Missing column: {c}")
    if df[id_col].isna().any() or df[id_col].duplicated().any():
        raise ValueError(f"{id_col} must be unique and non-missing.")
    markers=[c for c in df.columns if c not in {id_col,population_col}]
    if not markers: raise ValueError("No marker columns found.")
    try: x=df[markers].to_numpy(dtype=float)
    except Exception as e: raise ValueError("Marker values must be numeric 0/1/2/NaN.") from e
    obs=x[np.isfinite(x)]
    if not np.isin(obs,[0,1,2]).all():
        raise ValueError("Only diploid biallelic dosage 0, 1, 2 and NaN are supported.")
    return markers

def pop_groups(df,population_col="population_id"):
    if population_col not in df: raise ValueError(f"Missing {population_col}.")
    pops=sorted(df[population_col].dropna().astype(str).unique())
    if not pops: raise ValueError("No population labels found.")
    return [(p,df.loc[df[population_col].astype(str)==p]) for p in pops]

def marker_af(x):
    x=np.asarray(x,dtype=float); x=x[np.isfinite(x)]
    if len(x)==0: return np.nan,0
    return float(x.sum()/(2*len(x))),len(x)

def save(df,path):
    from pathlib import Path
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True); df.to_csv(p,index=False); return p

def allele_freqs_for_pop(df,markers,population_col):
    result={}
    for pop,sub in pop_groups(df,population_col):
        result[pop]=np.array([marker_af(sub[m])[0] for m in markers],dtype=float)
    return result
