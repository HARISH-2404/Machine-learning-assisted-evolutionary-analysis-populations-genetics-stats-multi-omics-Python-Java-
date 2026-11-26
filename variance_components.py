"""Balanced random-effects variance component estimators by method of moments."""
import math

def one_way_random_effects(groups):
    """Estimate genotype and residual variance from balanced replicate groups.

    groups: mapping genotype -> list of replicate observations.
    Requires equal replication and at least two genotypes and replicates.
    """
    keys=list(groups)
    if len(keys)<2: raise ValueError("At least two groups are required")
    vals=[list(map(float,groups[k])) for k in keys]
    r=len(vals[0])
    if r<2 or any(len(v)!=r for v in vals): raise ValueError("Balanced groups with >=2 replicates required")
    means=[sum(v)/r for v in vals]
    grand=sum(sum(v) for v in vals)/(len(vals)*r)
    ss_between=r*sum((m-grand)**2 for m in means)
    ss_within=sum(sum((x-m)**2 for x in v) for v,m in zip(vals,means))
    ms_between=ss_between/(len(vals)-1)
    ms_within=ss_within/(len(vals)*(r-1))
    vg=max(0.0,(ms_between-ms_within)/r)
    return {"genetic_variance":vg,"residual_variance":ms_within,
            "ms_genotype":ms_between,"ms_error":ms_within,
            "genotypes":len(vals),"replicates_per_genotype":r}

def plot_level_components(genetic_variance, genotype_by_environment_variance,
                          residual_variance, environments=1, replications=1):
    """Return variance components and plot-level phenotypic variance.

    Inputs are assumed already estimated for a compatible balanced model.
    """
    if environments<1 or replications<1: raise ValueError("Counts must be positive")
    vp=genetic_variance+genotype_by_environment_variance/environments+residual_variance/(environments*replications)
    return {"Vg":genetic_variance,"Vge":genotype_by_environment_variance,
            "Ve":residual_variance,"Vp_mean":vp}
