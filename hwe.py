"""Hardy-Weinberg equilibrium Pearson chi-square diagnostic (1 df)."""
import numpy as np
import pandas as pd
from math import erfc,sqrt
from ._core import validate_dosage,pop_groups,save

def hardy_weinberg(genotypes,population_col=None,id_col="genotype_id",
                   min_expected=5.0,output=None):
    markers=validate_dosage(genotypes,id_col,population_col)
    groups=[("ALL",genotypes)] if population_col is None else pop_groups(genotypes,population_col)
    rows=[]
    for pop,sub in groups:
        for m in markers:
            x=sub[m].to_numpy(float); x=x[np.isfinite(x)]; n=len(x)
            if n==0: continue
            p=x.sum()/(2*n); q=1-p
            obs=np.array([(x==0).sum(),(x==1).sum(),(x==2).sum()],float)
            exp=n*np.array([q*q,2*p*q,p*p])
            chi=float(np.sum((obs-exp)**2/exp)) if np.all(exp>0) else np.nan
            pv=erfc(sqrt(chi/2)) if np.isfinite(chi) else np.nan
            rows.append({"population_id":pop,"marker_id":m,"n":n,"p_alt":p,
                         "obs_00":int(obs[0]),"obs_01":int(obs[1]),"obs_11":int(obs[2]),
                         "exp_00":exp[0],"exp_01":exp[1],"exp_11":exp[2],
                         "chi_square":chi,"p_value_asymptotic":pv,
                         "expected_counts_all_ge_min":bool(np.all(exp>=min_expected))})
    out=pd.DataFrame(rows)
    if output: save(out,output)
    return out
