"""Shared transcriptomics matrix and metadata validation."""
import numpy as np
import pandas as pd

def count_matrix(data, require_integer=True):
    """Validate genes x samples count matrix and return float DataFrame."""
    df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
    if df.empty or df.shape[0]<1 or df.shape[1]<1: raise ValueError("count matrix must not be empty")
    x=df.to_numpy(dtype=float)
    if not np.all(np.isfinite(x)) or np.any(x<0): raise ValueError("counts must be finite and non-negative")
    if require_integer and not np.allclose(x,np.round(x)): raise ValueError("raw counts must be integers; TPM/FPKM are not raw counts")
    if df.index.has_duplicates or df.columns.has_duplicates: raise ValueError("gene and sample identifiers must be unique")
    return df.astype(float)

def validate_metadata(metadata, sample_ids):
    """Require one metadata row per sample and exact sample-ID matching."""
    md=metadata.copy() if isinstance(metadata,pd.DataFrame) else pd.DataFrame(metadata)
    if "sample_id" in md.columns: md=md.set_index("sample_id")
    if md.index.has_duplicates: raise ValueError("metadata sample IDs must be unique")
    missing=set(sample_ids)-set(md.index)
    if missing: raise ValueError(f"metadata missing samples: {sorted(missing)}")
    return md.loc[list(sample_ids)].copy()
