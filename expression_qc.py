"""Quality-control summaries for RNA-seq raw count matrices."""
import numpy as np
import pandas as pd
from ._core import count_matrix

def expression_qc(counts):
    df=count_matrix(counts)
    lib=df.sum(axis=0)
    detected=(df>0).sum(axis=0)
    zero_fraction=(df==0).mean(axis=0)
    rows=[]
    for s in df.columns:
        rows.append({"sample_id":s,"library_size":float(lib[s]),"detected_genes":int(detected[s]),
                     "zero_fraction":float(zero_fraction[s]),"median_gene_count":float(df[s].median())})
    return pd.DataFrame(rows)

def gene_qc(counts):
    df=count_matrix(counts)
    return pd.DataFrame({"total_count":df.sum(axis=1),"samples_detected":(df>0).sum(axis=1),
                         "mean_count":df.mean(axis=1),"variance":df.var(axis=1)})

def flag_library_size_outliers(counts, z_threshold=3.0):
    """Flag samples by absolute z-score of log library size; does not remove samples."""
    lib=count_matrix(counts).sum(axis=0)
    x=np.log1p(lib.to_numpy()); sd=x.std(ddof=0)
    z=np.zeros_like(x) if sd==0 else (x-x.mean())/sd
    return pd.DataFrame({"sample_id":lib.index,"library_size":lib.values,"log_library_z":z,
                         "outlier":np.abs(z)>z_threshold})
