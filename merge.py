"""Safe merging of matched omics feature blocks."""
import pandas as pd

def merge_omics_blocks(blocks, join="inner", prefixes=True):
    """Concatenate sample-by-feature blocks. Rows must be samples, columns features."""
    if not blocks: raise ValueError("No omics blocks supplied")
    if join not in {"inner", "outer"}: raise ValueError("join must be inner or outer")
    prepared=[]
    for name, df in blocks.items():
        if not df.index.is_unique: raise ValueError(f"Duplicate sample IDs in {name}")
        if not df.columns.is_unique: raise ValueError(f"Duplicate feature IDs in {name}")
        x=df.copy()
        x.index=x.index.astype(str)
        if prefixes: x.columns=[f"{name}::{c}" for c in x.columns]
        prepared.append(x)
    return pd.concat(prepared, axis=1, join=join, sort=False)
