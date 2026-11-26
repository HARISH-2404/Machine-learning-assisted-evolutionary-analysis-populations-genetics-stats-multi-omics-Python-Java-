"""ANOVA-based variance partitioning of diploid SNP dosage by population.

This is an explicitly labelled one-level, dosage-based AMOVA approximation.
For each marker, one-way random-effects ANOVA partitions within- and
among-population variance. For unequal sample sizes:
n0 = (N - sum(n_i^2)/N)/(k-1)
sigma2_between = (MS_between - MS_within)/n0
sigma2_within = MS_within
The among-population variance fraction is sigma2_between /
(sigma2_between+sigma2_within). This is not a substitute for a full
molecular-distance AMOVA (e.g. Arlequin/poppr) when hierarchical design,
haplotypes or molecular distances are required.
"""
import numpy as np
import pandas as pd
from ._core import validate_dosage,pop_groups,save

def amova(genotypes,population_col="population_id",id_col="genotype_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=dict(pop_groups(genotypes,population_col)); pops=list(groups)
    rows=[]
    for m in markers:
        vals=[]; labels=[]
        for pop,sub in groups.items():
            x=sub[m].dropna().to_numpy(float)
            vals.extend(x.tolist()); labels.extend([pop]*len(x))
        y=np.asarray(vals,float); lab=np.asarray(labels); N=len(y); k=len(set(lab))
        if N<=k or k<2: continue
        grand=y.mean()
        ssb=sum(np.sum(lab==p)*(y[lab==p].mean()-grand)**2 for p in set(lab))
        ssw=sum(np.sum((y[lab==p]-y[lab==p].mean())**2) for p in set(lab))
        dfb=k-1; dfw=N-k; msb=ssb/dfb; msw=ssw/dfw
        ns=np.array([np.sum(lab==p) for p in set(lab)],float)
        n0=(N-np.sum(ns**2)/N)/dfb
        vb=(msb-msw)/n0 if n0>0 else np.nan
        vw=msw
        # Negative component estimates can occur due to sampling noise.
        phi=vb/(vb+vw) if np.isfinite(vb) and vb+vw>0 else np.nan
        rows.append({"marker_id":m,"n":N,"populations":k,"SS_among":ssb,
                     "SS_within":ssw,"df_among":dfb,"df_within":dfw,
                     "MS_among":msb,"MS_within":msw,"variance_among":vb,
                     "variance_within":vw,"Phi_among":phi})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
