"""Allele frequency, minor-allele frequency and call-rate estimation."""
import numpy as np
import pandas as pd
from ._core import validate_dosage,pop_groups,marker_af,save

def allele_frequency(genotypes,population_col=None,id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=[("ALL",genotypes)] if population_col is None else pop_groups(genotypes,population_col)
    rows=[]
    for pop,sub in groups:
        for m in markers:
            x=sub[m].to_numpy(float); p,n=marker_af(x); q=1-p if np.isfinite(p) else np.nan
            rows.append({"population_id":pop,"marker_id":m,"n_called":n,
                         "call_rate":n/len(sub),"ref_frequency":q,"alt_frequency":p,
                         "maf":min(p,q) if np.isfinite(p) else np.nan})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
