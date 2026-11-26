"""RNA-seq preprocessing, normalization and sample QC helpers."""
import numpy as np
import pandas as pd
from ._core import count_matrix
from .normalization import cpm, median_ratio_size_factors

def summarize_rnaseq(counts):
    """Return library size, detected genes, zero fraction and count quantiles per sample."""
    df=count_matrix(counts)
    rows=[]
    for s in df.columns:
        x=df[s].to_numpy(float)
        rows.append({"sample_id":s,"library_size":float(x.sum()),
                     "detected_genes":int((x>0).sum()),"zero_fraction":float(np.mean(x==0)),
                     "mean_count":float(np.mean(x)),"median_count":float(np.median(x)),
                     "q90_count":float(np.quantile(x,.9))})
    return pd.DataFrame(rows)

def rna_seq_normalize(counts, method="cpm"):
    """Normalize counts using CPM or DESeq-style median-ratio size factors."""
    if method=="cpm": return cpm(counts), None
    if method in ("median_ratio","deseq_size_factor"):
        normalized,sf=__import__("popgenml.transcriptomics.normalization",fromlist=["normalize_median_ratio"]).normalize_median_ratio(counts)
        return normalized,sf
    raise ValueError("method must be 'cpm' or 'median_ratio'")

def filter_expression(counts, min_count=10, min_samples=2):
    """Keep genes meeting a minimum count in a minimum number of samples."""
    df=count_matrix(counts)
    keep=(df>=min_count).sum(axis=1)>=min_samples
    return df.loc[keep].copy()
