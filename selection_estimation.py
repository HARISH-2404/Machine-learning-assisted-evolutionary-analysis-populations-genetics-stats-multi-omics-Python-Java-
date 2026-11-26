"""Approximate selection coefficient estimation from temporal allele frequencies."""
import numpy as np
import pandas as pd

def estimate_selection_logit(trajectory, group_cols=("chromosome","position","gene_id","allele","population","location"),
                             generation_col="generation", frequency_col="allele_frequency",
                             epsilon=1e-4):
    """Estimate a descriptive logit-frequency slope per generation.

    Under a simple deterministic haploid/additive approximation, logit(p_t)
    changes approximately linearly with generation, with slope related to s.
    This is NOT a general diploid selection estimator. Results are sensitive to
    drift, migration, structure, sampling error and dominance.
    """
    rows=[]
    for key, d in trajectory.groupby(list(group_cols),dropna=False):
        d=d.sort_values(generation_col)
        p=np.clip(d[frequency_col].to_numpy(float),epsilon,1-epsilon)
        gen=d[generation_col].to_numpy(float)
        if len(np.unique(gen))<2:
            slope=np.nan; r2=np.nan
        else:
            y=np.log(p/(1-p))
            slope,intercept=np.polyfit(gen,y,1)
            pred=intercept+slope*gen
            ssr=((y-pred)**2).sum(); sst=((y-y.mean())**2).sum()
            r2=1-ssr/sst if sst>0 else np.nan
        key=key if isinstance(key,tuple) else (key,)
        rows.append(dict(zip(group_cols,key)) | {"selection_slope_approx":slope,
                    "trajectory_r_squared":r2,"n_timepoints":len(d)})
    return pd.DataFrame(rows)

def deterministic_selection_step(p, s, dominance=0.5):
    """One-generation diploid biallelic viability-selection update.

    Fitness: AA=1+s, Aa=1+h*s, aa=1. Requires p=frequency of A.
    """
    if not 0<=p<=1: raise ValueError("p must be in [0,1]")
    q=1-p
    wAA=1+s; wAa=1+dominance*s; waa=1.0
    if min(wAA,wAa,waa)<0: raise ValueError("Fitness values cannot be negative")
    mean_w=p*p*wAA+2*p*q*wAa+q*q*waa
    if mean_w==0: raise ValueError("Mean fitness is zero")
    return (p*p*wAA+p*q*wAa)/mean_w

def simulate_selection_trajectory(p0, generations, s, dominance=0.5):
    if generations<0: raise ValueError("generations must be nonnegative")
    out=[float(p0)]
    for _ in range(generations): out.append(deterministic_selection_step(out[-1],s,dominance))
    return pd.DataFrame({"generation":range(generations+1),"allele_frequency":out})
