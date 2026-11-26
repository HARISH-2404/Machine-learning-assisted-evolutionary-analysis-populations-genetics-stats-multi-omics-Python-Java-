"""Read and summarize gene-level count matrices."""
import pandas as pd
from ._core import count_matrix

def read_count_matrix(path, sep="\t", gene_id_column=0):
    df=pd.read_csv(path,sep=sep,index_col=gene_id_column)
    return count_matrix(df,require_integer=True)

def library_sizes(counts):
    df=count_matrix(counts)
    return df.sum(axis=0).rename("library_size").to_frame()

def filter_low_count_genes(counts, min_count=10, min_samples=2):
    """Keep genes with >= min_count in at least min_samples samples."""
    df=count_matrix(counts)
    if min_count<0 or min_samples<1: raise ValueError("invalid filtering thresholds")
    keep=(df>=min_count).sum(axis=1)>=min_samples
    return df.loc[keep].copy()
