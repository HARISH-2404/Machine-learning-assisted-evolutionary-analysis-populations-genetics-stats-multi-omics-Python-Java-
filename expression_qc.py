"""RNA-seq count-matrix preprocessing with library sizes, CPM and log2(CPM+1).

Requires genuine nonnegative integer raw counts for use as an RNA-seq count
preprocessing step. Do not feed TPM, FPKM, log-transformed expression, or
continuous synthetic expression into a count-based differential-expression
model. CPM is for exploration, not a substitute for DESeq2 size factors.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._utils import require_unique,numeric_matrix,save_table

class ExpressionQC:
    def __init__(self,min_total_count=10,min_samples_expressed=3,
                 min_library_size=1):
        if min_total_count<0 or min_samples_expressed<1 or min_library_size<0:
            raise ValueError("Count and sample thresholds must be nonnegative.")
        self.min_total_count=min_total_count
        self.min_samples_expressed=min_samples_expressed
        self.min_library_size=min_library_size

    def process(self,frame,sample_col="sample_id",metadata_cols=("group",),
                output_dir=None):
        require_unique(frame,sample_col)
        metadata=[c for c in metadata_cols if c in frame.columns]
        genes=[c for c in frame.columns if c not in [sample_col]+metadata]
        if not genes: raise ValueError("No gene columns found.")
        x=numeric_matrix(frame,genes)
        if not np.isfinite(x).all() or (x<0).any() or not np.all(x==np.floor(x)):
            raise ValueError("Expression counts must be finite, nonnegative integers.")
        lib=x.sum(axis=1)
        sample_qc=pd.DataFrame({sample_col:frame[sample_col].to_numpy(),
                                "library_size":lib,"detected_genes":(x>0).sum(axis=1),
                                "retained":lib>=self.min_library_size})
        keep_samples=sample_qc.retained.to_numpy()
        if not keep_samples.any(): raise ValueError("No samples passed library-size QC.")
        x2=x[keep_samples]
        sample_frame=frame.loc[keep_samples,[sample_col]+metadata].reset_index(drop=True)
        gene_keep=(x2.sum(axis=0)>=self.min_total_count)&(
            (x2>0).sum(axis=0)>=self.min_samples_expressed)
        gene_qc=pd.DataFrame({"gene_id":genes,"total_count":x2.sum(axis=0),
                              "samples_expressed":(x2>0).sum(axis=0),
                              "retained":gene_keep})
        if not gene_keep.any(): raise ValueError("No genes passed expression QC.")
        x3=x2[:,gene_keep].astype(np.int64)
        filtered=pd.concat([sample_frame,pd.DataFrame(x3,columns=np.array(genes)[gene_keep])],axis=1)
        # CPM denominator uses total library counts BEFORE gene filtering.
        denominators=lib[keep_samples]
        cpm=x3/np.maximum(denominators[:,None],1)*1_000_000
        normalized=pd.concat([sample_frame,pd.DataFrame(np.log2(cpm+1),
                             columns=np.array(genes)[gene_keep])],axis=1)
        if output_dir:
            save_table(filtered,f"{output_dir}/expression_filtered_raw_counts.csv")
            save_table(normalized,f"{output_dir}/expression_log2_cpm_plus1.csv")
            save_table(sample_qc,f"{output_dir}/expression_sample_qc.csv")
            save_table(gene_qc,f"{output_dir}/expression_gene_qc.csv")
        return filtered,normalized,sample_qc,gene_qc
