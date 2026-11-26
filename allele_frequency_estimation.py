"""Allele-frequency, heterozygosity and temporal change estimators."""
import numpy as np
import pandas as pd

def estimate_allele_frequency(counts, allele_count_col="allele_count", total_col="total_alleles"):
    x=counts.copy()
    if (x[total_col] <= 0).any() or (x[allele_count_col] < 0).any() or (x[allele_count_col] > x[total_col]).any():
        raise ValueError("Require 0 <= allele_count <= total_alleles and total_alleles > 0")
    x["allele_frequency"]=x[allele_count_col]/x[total_col]
    x["other_allele_frequency"]=1-x["allele_frequency"]
    x["expected_heterozygosity"]=2*x["allele_frequency"]*x["other_allele_frequency"]
    return x

def temporal_frequency_change(data, group_cols=("chromosome","position","gene_id","allele","population","location"),
                              generation_col="generation", frequency_col="allele_frequency"):
    x=data.sort_values([*group_cols,generation_col]).copy()
    g=x.groupby(list(group_cols),dropna=False)
    x["previous_frequency"]=g[frequency_col].shift(1)
    x["previous_generation"]=g[generation_col].shift(1)
    x["delta_frequency"]=x[frequency_col]-x["previous_frequency"]
    delta_gen=x[generation_col]-x["previous_generation"]
    x["delta_frequency_per_generation"]=x["delta_frequency"]/delta_gen.replace(0,np.nan)
    return x

def estimate_fis(observed_heterozygosity, allele_frequency):
    p=np.asarray(allele_frequency,dtype=float)
    ho=np.asarray(observed_heterozygosity,dtype=float)
    he=2*p*(1-p)
    return np.where(he>0,1-ho/he,np.nan)
