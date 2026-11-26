"""Observed and projected quantitative phenotype trajectories."""
import pandas as pd

def summarize_phenotype_trajectory(data, trait_cols, generation_col="generation",
                                   population_col="population"):
    required={generation_col,population_col,*trait_cols}
    if not required.issubset(data.columns):
        raise ValueError(f"Missing columns: {sorted(required-set(data.columns))}")
    return (data.groupby([population_col,generation_col])[trait_cols]
            .agg(["mean","std","count"]).reset_index())

def project_phenotype_trajectory(allele_paths, baseline_mean, additive_effect,
                                 baseline_frequency=None, dominance_effect=0.0):
    """Project trait mean from one locus under HW and supplied additive/dominance effects."""
    x=allele_paths.copy()
    p=x["allele_frequency"].astype(float)
    p0=float(p.iloc[0]) if baseline_frequency is None else float(baseline_frequency)
    if not 0<=p0<=1: raise ValueError("baseline_frequency must be in [0,1]")
    genetic=2*additive_effect*p+2*dominance_effect*p*(1-p)
    base=2*additive_effect*p0+2*dominance_effect*p0*(1-p0)
    x["projected_phenotype_mean"]=baseline_mean+genetic-base
    return x
