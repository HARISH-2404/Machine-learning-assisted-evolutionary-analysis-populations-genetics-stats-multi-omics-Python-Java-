"""Genotype QC for diploid biallelic SNP dosage (0, 1, 2).

Missing calls are NaN; quality filters include marker/sample call rate,
minor-allele frequency (MAF), heterozygosity and Hardy-Weinberg equilibrium.
The HWE test is a chi-square approximation for a single randomly mating
diploid population; avoid interpreting it in structured/inbred cotton panels.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._utils import save_table,require_unique,numeric_matrix

class GenotypeQC:
    def __init__(self,min_marker_call_rate=0.90,min_sample_call_rate=0.90,
                 min_maf=0.01,max_heterozygosity=1.0):
        for name,value in [("min_marker_call_rate",min_marker_call_rate),
                           ("min_sample_call_rate",min_sample_call_rate),
                           ("min_maf",min_maf),("max_heterozygosity",max_heterozygosity)]:
            if not 0<=value<=1:
                raise ValueError(f"{name} must be between 0 and 1.")
        if min_maf>0.5:
            raise ValueError("min_maf cannot exceed 0.5.")
        self.min_marker_call_rate=min_marker_call_rate
        self.min_sample_call_rate=min_sample_call_rate
        self.min_maf=min_maf
        self.max_heterozygosity=max_heterozygosity

    @staticmethod
    def marker_summary(frame,id_col="genotype_id"):
        require_unique(frame,id_col)
        markers=[c for c in frame if c!=id_col]
        if not markers: raise ValueError("No SNP columns found.")
        x=numeric_matrix(frame,markers)
        valid=np.isfinite(x)
        if not np.isin(x[valid],[0,1,2]).all():
            raise ValueError("Diploid SNP dosage must be 0, 1, 2 or missing.")
        called=valid.sum(axis=0)
        alt=np.divide(np.where(valid,x,0).sum(axis=0),2*called,
                      out=np.full(len(markers),np.nan),where=called>0)
        maf=np.minimum(alt,1-alt)
        hetero=np.divide(((x==1)&valid).sum(axis=0),called,
                         out=np.full(len(markers),np.nan),where=called>0)
        obs0=((x==0)&valid).sum(axis=0)
        obs1=((x==1)&valid).sum(axis=0)
        obs2=((x==2)&valid).sum(axis=0)
        # Pearson chi-square with 1 df, using the sample allele frequency.
        expected=np.stack([called*(1-alt)**2,called*2*alt*(1-alt),called*alt**2],axis=1)
        observed=np.stack([obs0,obs1,obs2],axis=1)
        chi2=np.divide((observed-expected)**2,expected,
                       out=np.zeros_like(expected,dtype=float),where=expected>0).sum(axis=1)
        # Chi-square(1) upper-tail probability = erfc(sqrt(chi2/2)).
        from math import erfc,sqrt
        hwe_p=np.array([erfc(sqrt(v/2)) if n>=3 and np.isfinite(v) and 0<p<1
                        else np.nan for v,n,p in zip(chi2,called,alt)])
        return pd.DataFrame({"marker_id":markers,"call_rate":called/len(frame),
                             "alt_allele_frequency":alt,"maf":maf,
                             "heterozygosity":hetero,"hwe_chi2_approx_p":hwe_p,
                             "n_called":called})

    @staticmethod
    def sample_summary(frame,id_col="genotype_id"):
        require_unique(frame,id_col)
        markers=[c for c in frame if c!=id_col]
        x=numeric_matrix(frame,markers)
        valid=np.isfinite(x)
        return pd.DataFrame({id_col:frame[id_col].to_numpy(),
                             "call_rate":valid.mean(axis=1),
                             "heterozygosity":np.divide(((x==1)&valid).sum(axis=1),
                                 valid.sum(axis=1),out=np.full(len(frame),np.nan),
                                 where=valid.sum(axis=1)>0)})

    def filter(self,frame,id_col="genotype_id",impute=True,output_dir=None):
        """Apply sample QC then marker QC; optionally impute marker modes.

        Returns filtered dosage, marker QC, sample QC, and a summary dictionary.
        HWE is reported for diagnostics but not automatically filtered.
        """
        sample_qc=self.sample_summary(frame,id_col)
        sample_keep=sample_qc.call_rate>=self.min_sample_call_rate
        sample_qc["retained"]=sample_keep.to_numpy()
        sub=frame.loc[sample_keep.to_numpy()].copy()
        if sub.empty: raise ValueError("No samples passed sample call-rate QC.")
        marker_qc=self.marker_summary(sub,id_col)
        keep=(marker_qc.call_rate>=self.min_marker_call_rate)&(
            marker_qc.maf>=self.min_maf)&(
            marker_qc.heterozygosity<=self.max_heterozygosity)
        marker_qc["retained"]=keep.to_numpy()
        kept=marker_qc.loc[keep,"marker_id"].tolist()
        if not kept: raise ValueError("No markers passed marker QC.")
        filtered=sub[[id_col]+kept].copy()
        if impute:
            for marker in kept:
                values=filtered[marker].dropna()
                if values.empty: raise ValueError(f"Cannot impute entirely missing marker {marker}.")
                # Most frequent observed genotype; deterministic tie-break to smaller dosage.
                mode=values.value_counts().sort_index().sort_values(ascending=False,kind="stable").index[0]
                filtered[marker]=filtered[marker].fillna(mode).astype(int)
        summary={"input_samples":len(frame),"retained_samples":len(sub),
                 "input_markers":len(frame.columns)-1,"retained_markers":len(kept),
                 "imputed":bool(impute)}
        if output_dir:
            save_table(filtered,f"{output_dir}/genotype_qc_filtered.csv")
            save_table(marker_qc,f"{output_dir}/marker_qc_summary.csv")
            save_table(sample_qc,f"{output_dir}/sample_qc_summary.csv")
        return filtered,marker_qc,sample_qc,summary
