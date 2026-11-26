"""Observed and expected heterozygosity per marker/population."""
import numpy as np
import pandas as pd
from ._core import validate_dosage,pop_groups,marker_af,save

def heterozygosity(genotypes,population_col=None,id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=[("ALL",genotypes)] if population_col is None else pop_groups(genotypes,population_col)
    rows=[]
    for pop,sub in groups:
        for m in markers:
            x=sub[m].to_numpy(float); x=x[np.isfinite(x)]; n=len(x)
            if n==0: ho=he=he_unbiased=np.nan; p=np.nan
            else:
                p=x.sum()/(2*n); q=1-p
                ho=float(np.mean(x==1)); he=2*p*q
                he_unbiased=he*n/(n-1) if n>1 else np.nan
            rows.append({"population_id":pop,"marker_id":m,"n_called":n,
                         "observed_heterozygosity":ho,"expected_heterozygosity":he,
                         "unbiased_expected_heterozygosity":he_unbiased})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
