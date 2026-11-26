"""Library-size normalization helpers for exploratory visualization.

These methods do not replace DESeq2 median-ratio size factors or edgeR TMM
for formal differential-expression inference.
"""
import numpy as np
import pandas as pd
from ._core import count_matrix

def cpm(counts, log=False, prior_count=0.5):
    df=count_matrix(counts); lib=df.sum(axis=0)
    if (lib<=0).any(): raise ValueError("zero-size library cannot be normalized")
    values=df.div(lib,axis=1)*1e6
    return np.log2(values+prior_count) if log else values

def median_ratio_size_factors(counts):
    """DESeq-style median-of-ratios size factors, excluding zero geometric means."""
    df=count_matrix(counts).to_numpy(dtype=float)
    with np.errstate(divide="ignore",invalid="ignore"):
        gm=np.exp(np.mean(np.log(np.where(df>0,df,np.nan)),axis=1))
    valid=np.isfinite(gm)&(gm>0)
    if not valid.any(): raise ValueError("no genes with positive counts in every sample")
    ratios=df[valid,:]/gm[valid,None]
    sf=np.nanmedian(ratios,axis=0)
    if np.any(~np.isfinite(sf)) or np.any(sf<=0): raise ValueError("could not estimate positive size factors")
    sf=sf/np.exp(np.mean(np.log(sf)))
    return pd.Series(sf,index=None,name="size_factor")

def normalize_median_ratio(counts):
    df=count_matrix(counts); sf=median_ratio_size_factors(df)
    sf.index=df.columns
    return df.div(sf,axis=1),sf
