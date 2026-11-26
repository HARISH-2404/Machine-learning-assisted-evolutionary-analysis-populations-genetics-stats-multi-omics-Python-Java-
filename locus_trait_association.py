"""Locus-level genotype–trait association screens."""
import numpy as np
import pandas as pd
from scipy.stats import linregress

def scan_locus_trait(genotype_dosage, phenotype, min_n=8):
    """Univariate dosage regression for each marker; returns beta, p and R2.

    Genotype dosage should encode allele copies consistently. This is a
    screening association, not a mixed-model GWAS; account for kinship,
    population structure, environment and multiple testing in formal analyses.
    """
    y=phenotype.astype(float)
    rows=[]
    for marker in genotype_dosage.columns:
        d=pd.concat([genotype_dosage[marker],y],axis=1).dropna()
        if len(d)<min_n or d.iloc[:,0].nunique()<2: continue
        fit=linregress(d.iloc[:,0].astype(float),d.iloc[:,1])
        rows.append({"marker":marker,"n":len(d),"effect_per_allele":fit.slope,
                     "intercept":fit.intercept,"r_squared":fit.rvalue**2,"p_value":fit.pvalue})
    out=pd.DataFrame(rows)
    if len(out):
        order=out.p_value.rank(method="min")
        out["fdr_bh"]=np.minimum(1,out.p_value*len(out)/order)
    return out.sort_values("p_value") if len(out) else out
