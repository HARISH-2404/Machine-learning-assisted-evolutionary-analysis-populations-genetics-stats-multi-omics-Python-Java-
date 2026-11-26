"""Practical quantitative-genetic estimators for balanced trial data."""
import numpy as np
import pandas as pd

def estimate_broad_sense_h2(data, genotype_col, environment_col, trait_col):
    """Balanced-design ANOVA moment estimator: H2=Vg/(Vg+Ve/r).

    Assumes genotypes evaluated in each environment with replication. This
    simplified estimator is not a substitute for REML in unbalanced/multi-year
    trials or complex spatial designs.
    """
    means=data.groupby([genotype_col,environment_col])[trait_col].mean().unstack()
    grand=data[trait_col].mean()
    gmeans=data.groupby(genotype_col)[trait_col].mean()
    n_env=data[environment_col].nunique()
    r=data.groupby([genotype_col,environment_col]).size().median()
    ms_g=n_env*r*((gmeans-grand)**2).sum()/(len(gmeans)-1) if len(gmeans)>1 else np.nan
    residual=data.merge(means.stack().rename("_cellmean").reset_index(),on=[genotype_col,environment_col])
    sse=((residual[trait_col]-residual["_cellmean"])**2).sum()
    df_e=max(len(data)-means.size,1)
    ms_e=sse/df_e
    vg=max((ms_g-ms_e)/(n_env*r),0) if np.isfinite(ms_g) else np.nan
    vp=vg+ms_e/r
    return {"genotypic_variance":vg,"residual_variance":ms_e,"broad_sense_h2":vg/vp if vp>0 else np.nan,
            "n_genotypes":len(gmeans),"n_environments":n_env,"replication_approx":r}

def estimate_genetic_correlation(genotype_trait_means, trait1, trait2, genotype_col="genotype"):
    x=genotype_trait_means[[genotype_col,trait1,trait2]].dropna()
    if len(x)<3: return np.nan
    return float(x[trait1].corr(x[trait2]))

def selection_differential(data, trait_col, selected_col):
    selected=data.loc[data[selected_col].astype(bool),trait_col].mean()
    base=data[trait_col].mean()
    return {"population_mean":float(base),"selected_mean":float(selected),
            "selection_differential":float(selected-base)}

def realized_response(parent_mean, offspring_mean):
    return float(offspring_mean-parent_mean)

def response_per_generation(response, generations):
    if generations<=0: raise ValueError("generations must be positive")
    return float(response/generations)

def additive_dominance_from_genotype_means(means, genotype_col="genotype_class", trait_col="trait_mean"):
    """Estimate a and d from aa/Aa/AA means, where class labels are 0,1,2 copies A."""
    m=means.set_index(genotype_col)[trait_col]
    if not {0,1,2}.issubset(set(m.index)): raise ValueError("Need genotype classes 0,1,2")
    aa=float(m.loc[0]); het=float(m.loc[1]); AA=float(m.loc[2])
    return {"additive_effect_a":(AA-aa)/2,"dominance_deviation_d":het-(AA+aa)/2}
