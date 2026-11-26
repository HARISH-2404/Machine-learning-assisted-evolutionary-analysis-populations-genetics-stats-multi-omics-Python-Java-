"""Observed genotype frequencies for dosage classes 0, 1 and 2."""
import numpy as np
import pandas as pd
from ._core import validate_dosage,pop_groups,save

def genotype_frequency(genotypes,population_col=None,id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=[("ALL",genotypes)] if population_col is None else pop_groups(genotypes,population_col)
    rows=[]
    for pop,sub in groups:
        for m in markers:
            x=sub[m].to_numpy(float); x=x[np.isfinite(x)]; n=len(x)
            counts=[int((x==i).sum()) for i in (0,1,2)]
            rows.append({"population_id":pop,"marker_id":m,"n_called":n,
                         "n_0":counts[0],"n_1":counts[1],"n_2":counts[2],
                         "freq_0":counts[0]/n if n else np.nan,
                         "freq_1":counts[1]/n if n else np.nan,
                         "freq_2":counts[2]/n if n else np.nan})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
