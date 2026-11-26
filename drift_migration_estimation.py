"""Temporal drift, effective-size and migration estimators."""
import numpy as np
import pandas as pd

def estimate_temporal_ne(trajectory, group_cols=("chromosome","position","gene_id","allele","population","location"),
                         generation_col="generation", frequency_col="allele_frequency"):
    """Rough diploid temporal Ne estimator using allele-frequency variance.

    For a neutral Wright-Fisher approximation over t generations:
    Var(delta p) approx p0(1-p0)*t/(2Ne). This estimator is noisy and should
    aggregate independent neutral loci; do not infer Ne from one selected locus.
    """
    rows=[]
    for key,d in trajectory.groupby(list(group_cols),dropna=False):
        d=d.sort_values(generation_col)
        if len(d)<2: continue
        p0=float(d[frequency_col].iloc[0]); p1=float(d[frequency_col].iloc[-1])
        t=float(d[generation_col].iloc[-1]-d[generation_col].iloc[0])
        denom=2*(p1-p0)**2
        ne=(p0*(1-p0)*t/denom) if denom>0 and t>0 else np.nan
        key=key if isinstance(key,tuple) else (key,)
        rows.append(dict(zip(group_cols,key)) | {"temporal_ne_rough":ne,"generations_elapsed":t})
    return pd.DataFrame(rows)

def estimate_migration_from_admixture(p_local, p_source, p_observed):
    """Single-source mixture estimate m=(p_observed-p_local)/(p_source-p_local)."""
    if abs(p_source-p_local)<1e-12: return np.nan
    return float(np.clip((p_observed-p_local)/(p_source-p_local),0,1))

def heterozygosity_decay(heterozygosity, generations, ne):
    """Neutral expected H_t = H_0*(1-1/(2Ne))^t for diploids."""
    if ne<=0 or generations<0: raise ValueError("Ne > 0 and generations >= 0 required")
    return float(heterozygosity*(1-1/(2*ne))**generations)
