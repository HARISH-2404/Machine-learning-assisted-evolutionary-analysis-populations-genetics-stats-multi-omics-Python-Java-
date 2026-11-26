"""Polymorphism information content (PIC) for biallelic SNP markers.

For allele frequencies p and q=1-p, Botstein's biallelic PIC is
PIC = 1 - (p^2+q^2) - 2p^2q^2 = 2pq(1-pq).
"""
import numpy as np
import pandas as pd
from ._core import validate_dosage,pop_groups,save

def polymorphism_information_content(genotypes,population_col=None,
                                     id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=[("ALL",genotypes)] if population_col is None else pop_groups(genotypes,population_col)
    rows=[]
    for pop,sub in groups:
        for m in markers:
            x=sub[m].dropna().to_numpy(float); n=len(x)
            p=x.sum()/(2*n) if n else np.nan; q=1-p if n else np.nan
            pic=2*p*q*(1-p*q) if n else np.nan
            rows.append({"population_id":pop,"marker_id":m,"n_called":n,
                         "p_alt":p,"q_ref":q,"PIC":pic})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
