"""Forward projection of quantitative traits under changing allele frequency."""
import pandas as pd

def project_trait_from_allele_frequency(paths, baseline_mean, allele_effect,
                                        baseline_frequency=None, dominance_effect=0.0):
    """Project mean trait using additive + heterozygote-dominance contribution.

    Allele effect is per-copy additive effect a; dominance_effect is d.
    Assumes Hardy-Weinberg genotype proportions and one causal biallelic locus.
    This is a transparent single-locus projection, not a genome-wide breeding-value model.
    """
    x=paths.copy()
    p=x["allele_frequency"].astype(float); q=1-p
    if baseline_frequency is None: baseline_frequency=float(p.iloc[0])
    p0=baseline_frequency; q0=1-p0
    genetic_mean=2*allele_effect*p + 2*dominance_effect*p*q
    base_genetic=2*allele_effect*p0+2*dominance_effect*p0*q0
    x["projected_trait_mean"]=baseline_mean+(genetic_mean-base_genetic)
    return x

def project_polygenic_trait(paths, locus_effects, baseline_mean, baseline_frequencies):
    """Sum independent locus contributions; locus_effects maps locus to additive effect."""
    x=paths.copy()
    total=0.0
    for locus,effect in locus_effects.items():
        col=f"p::{locus}"
        if col not in x: raise ValueError(f"Missing allele-frequency column {col}")
        p=x[col].astype(float); p0=float(baseline_frequencies[locus])
        total=total+2*effect*(p-p0)
    x["projected_polygenic_trait_mean"]=baseline_mean+total
    return x
