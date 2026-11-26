"""Nei's standard genetic distance from allele-frequency vectors.

For each population pair, I = sum_l sum_a sqrt(p_ial p_jal) /
sqrt[(sum_l sum_a p_ial^2)(sum_l sum_a p_jal^2)].
For biallelic loci, both allele frequencies are included. D = -ln(I).
Only loci with finite frequencies in both populations are used.
"""
import numpy as np
import pandas as pd
from itertools import combinations
from ._core import validate_dosage,pop_groups,save

def nei_genetic_distance(genotypes,population_col="population_id",
                         id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=dict(pop_groups(genotypes,population_col)); pops=list(groups)
    freqs={}
    for pop,sub in groups.items():
        p=[]
        for m in markers:
            x=sub[m].dropna().to_numpy(float)
            p.append(x.sum()/(2*len(x)) if len(x) else np.nan)
        freqs[pop]=np.asarray(p)
    rows=[]
    for a,b in combinations(pops,2):
        p1,p2=freqs[a],freqs[b]
        valid=np.isfinite(p1)&np.isfinite(p2)
        if not valid.any(): ident=dist=np.nan
        else:
            v1=np.column_stack([p1[valid],1-p1[valid]]).ravel()
            v2=np.column_stack([p2[valid],1-p2[valid]]).ravel()
            den=np.sqrt(np.sum(v1*v1)*np.sum(v2*v2))
            ident=float(np.sum(v1*v2)/den) if den>0 else np.nan
            dist=float(-np.log(max(ident,1e-15))) if np.isfinite(ident) else np.nan
        rows.append({"population_1":a,"population_2":b,"nei_identity":ident,
                     "nei_standard_distance":dist,"loci_used":int(valid.sum())})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
