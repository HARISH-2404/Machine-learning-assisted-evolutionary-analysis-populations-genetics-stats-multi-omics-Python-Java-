"""Early- and late-integration utilities."""
import numpy as np
import pandas as pd

def early_integration_matrix(blocks, sample_order=None, join="inner"):
    """Return sample-by-feature concatenated matrix and aligned sample IDs."""
    from .merge import merge_omics_blocks
    merged=merge_omics_blocks(blocks,join=join,prefixes=True)
    if sample_order is not None: merged=merged.reindex(sample_order)
    return merged

def late_integration_scores(score_tables, weights=None):
    """Weighted mean of per-omics prediction scores aligned by sample index."""
    if not score_tables: raise ValueError("No score tables")
    names=list(score_tables)
    w=np.ones(len(names),dtype=float) if weights is None else np.array([weights[n] for n in names],dtype=float)
    if np.any(w<0) or w.sum()==0: raise ValueError("Weights must be nonnegative with positive sum")
    aligned=pd.concat([score_tables[n].rename(n) for n in names],axis=1,join="inner")
    return pd.Series(aligned.to_numpy() @ (w/w.sum()),index=aligned.index,name="integrated_score")
