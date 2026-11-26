"""Transcript isoform abundance and usage calculations."""
import numpy as np
import pandas as pd

def isoform_usage(abundance, transcript_to_gene):
    """Convert transcript abundance (transcripts x samples) to within-gene proportions.

    Input may be estimated counts or TPM, but all isoforms in a gene must use
    the same abundance measure. Output is a long-form DataFrame.
    """
    df=abundance.copy() if isinstance(abundance,pd.DataFrame) else pd.DataFrame(abundance)
    if df.index.has_duplicates or np.any(df.to_numpy(dtype=float)<0): raise ValueError("abundance IDs must be unique and values non-negative")
    missing=set(df.index)-set(transcript_to_gene)
    if missing: raise ValueError(f"transcripts missing gene mapping: {sorted(missing)[:5]}")
    long=df.copy(); long["gene_id"]=[transcript_to_gene[x] for x in long.index]
    long["transcript_id"]=long.index
    melted=long.melt(id_vars=["gene_id","transcript_id"],var_name="sample_id",value_name="abundance")
    totals=melted.groupby(["gene_id","sample_id"])["abundance"].transform("sum")
    melted["isoform_usage"]=np.divide(melted.abundance,totals,out=np.zeros(len(melted),float),where=totals.to_numpy()>0)
    return melted

def dominant_isoform(usage_table):
    df=pd.DataFrame(usage_table)
    idx=df.groupby(["gene_id","sample_id"])["isoform_usage"].idxmax()
    return df.loc[idx,["gene_id","sample_id","transcript_id","isoform_usage"]].reset_index(drop=True)
