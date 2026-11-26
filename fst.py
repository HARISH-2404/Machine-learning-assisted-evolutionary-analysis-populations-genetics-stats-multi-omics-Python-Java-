"""Pairwise Hudson-style FST, ratio of summed locus numerators/denominators.

Per locus numerator=(p1-p2)^2 - p1(1-p1)/(n1-1)-p2(1-p2)/(n2-1).
Denominator=p1(1-p2)+p2(1-p1). Negative contributions are retained.
"""
import numpy as np
import pandas as pd
from itertools import combinations
from ._core import validate_dosage,pop_groups,save

def pairwise_fst(genotypes,population_col="population_id",id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=dict(pop_groups(genotypes,population_col)); pops=list(groups)
    rows=[]
    for a,b in combinations(pops,2):
        nums=[]; dens=[]
        for m in markers:
            x=groups[a][m].dropna().to_numpy(float); y=groups[b][m].dropna().to_numpy(float)
            n1,n2=len(x),len(y)
            if n1<2 or n2<2: continue
            p1=x.sum()/(2*n1); p2=y.sum()/(2*n2)
            den=p1*(1-p2)+p2*(1-p1)
            if den<=0: continue
            num=(p1-p2)**2-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1)
            nums.append(num); dens.append(den)
        fst=sum(nums)/sum(dens) if dens and sum(dens)>0 else np.nan
        rows.append({"population_1":a,"population_2":b,"FST_Hudson":fst,
                     "markers_used":len(dens)})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
