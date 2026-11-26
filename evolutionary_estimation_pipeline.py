"""Orchestrate parameter estimation and create a locus/population evolution table."""
import pandas as pd
from .allele_frequency_estimation import estimate_allele_frequency, temporal_frequency_change
from .selection_estimation import estimate_selection_logit
from .drift_migration_estimation import estimate_temporal_ne
from .population_differentiation import pairwise_fst_from_frequencies, summarize_diversity

def estimate_evolutionary_parameters(allele_counts):
    """Run descriptive estimators and return named result tables.

    Input must include chromosome, position, gene_id, allele, population,
    location, timepoint, generation, allele_count and total_alleles.
    """
    freq=estimate_allele_frequency(allele_counts)
    trajectory=temporal_frequency_change(freq)
    return {
        "allele_trajectory":trajectory,
        "selection_trend":estimate_selection_logit(trajectory),
        "temporal_ne_rough":estimate_temporal_ne(trajectory),
        "pairwise_fst_approx":pairwise_fst_from_frequencies(freq),
        "population_diversity":summarize_diversity(freq)
    }

def integrate_evolutionary_traits(trajectory, trait_data, expression_data=None):
    """Join quantitative/phenotypic and optional expression observations.

    Trait table must include population, location, timepoint, generation.
    Expression table must include same keys plus gene_id and expression.
    """
    keys=["population","location","timepoint","generation"]
    for c in keys:
        if c not in trait_data: raise ValueError(f"trait_data missing {c}")
    out=trajectory.merge(trait_data,on=keys,how="left",validate="many_to_many")
    if expression_data is not None:
        needed=set(keys+["gene_id"])
        if not needed.issubset(expression_data.columns):
            raise ValueError(f"expression_data requires {sorted(needed)}")
        out=out.merge(expression_data,on=keys+["gene_id"],how="left",validate="many_to_many")
    return out
