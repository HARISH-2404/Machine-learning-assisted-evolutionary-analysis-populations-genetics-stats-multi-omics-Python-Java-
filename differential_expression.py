"""Simple two-group differential-expression screening using log-CPM Welch tests.

For publication-quality RNA-seq DE, use DESeq2 or edgeR with biological
replicates and a design formula. This module is an educational baseline,
not a replacement for negative-binomial modelling.
"""
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from ._core import count_matrix, validate_metadata
from .normalization import cpm

def _bh(p):
    p=np.asarray(p,float); out=np.full(p.shape,np.nan)
    idx=np.where(np.isfinite(p)&(p>=0)&(p<=1))[0]
    if len(idx):
        order=idx[np.argsort(p[idx])]; m=len(order)
        q=p[order]*m/np.arange(1,m+1); q=np.minimum.accumulate(q[::-1])[::-1]
        out[order]=np.minimum(q,1)
    return out

def two_group_expression_test(counts, metadata, group_column, control, treatment,
                               min_total_count=10, log_prior=0.5):
    """Return log2FC, Welch t-test p and BH q for treatment versus control."""
    df=count_matrix(counts)
    md=validate_metadata(metadata,df.columns)
    if group_column not in md: raise ValueError(f"metadata lacks {group_column}")
    groups=md[group_column].astype(str)
    c=groups[groups==str(control)].index.tolist(); t=groups[groups==str(treatment)].index.tolist()
    if len(c)<2 or len(t)<2: raise ValueError("at least two biological samples per group are required")
    keep=df.sum(axis=1)>=min_total_count; df=df.loc[keep]
    norm=cpm(df,log=True,prior_count=log_prior)
    rows=[]
    for gene,row in norm.iterrows():
        a=row[t].to_numpy(); b=row[c].to_numpy()
        stat,p=ttest_ind(a,b,equal_var=False,nan_policy="omit")
        rows.append({"gene_id":gene,"mean_control":float(np.mean(b)),"mean_treatment":float(np.mean(a)),
                     "log2_fold_change":float(np.mean(a)-np.mean(b)),"t_statistic":float(stat),"p_value":float(p)})
    out=pd.DataFrame(rows)
    if len(out): out["p_adjusted_BH"]=_bh(out.p_value.to_numpy())
    return out.sort_values("p_value",na_position="last").reset_index(drop=True)
