"""Locus-aware allele-frequency trajectory analysis and forecasting."""
import pandas as pd
from .wright_fisher import simulate_wright_fisher, summarize_simulations

def build_allele_trajectory(observations, group_cols=("chromosome","position","gene_id","allele","population","location"),
                            generation_col="generation", count_col="allele_count", total_col="total_alleles"):
    required=set(group_cols)|{generation_col,count_col,total_col}
    if not required.issubset(observations.columns):
        raise ValueError(f"Missing columns: {sorted(required-set(observations.columns))}")
    x=observations.copy()
    if (x[total_col]<=0).any() or (x[count_col]<0).any() or (x[count_col]>x[total_col]).any():
        raise ValueError("Invalid allele counts")
    x["allele_frequency"]=x[count_col]/x[total_col]
    x=x.sort_values([*group_cols,generation_col])
    g=x.groupby(list(group_cols),dropna=False)
    x["previous_frequency"]=g["allele_frequency"].shift(1)
    x["frequency_change"]=x["allele_frequency"]-x["previous_frequency"]
    x["previous_generation"]=g[generation_col].shift(1)
    x["change_per_generation"]=x["frequency_change"]/(x[generation_col]-x["previous_generation"]).replace(0,pd.NA)
    return x

def forecast_locus_trajectory(p0, generations, population_size, selection_coefficient,
                              replicates=500, dominance=0.5, seed=1, **kwargs):
    paths=simulate_wright_fisher(p0,generations,population_size,
        selection_coefficient=selection_coefficient,dominance=dominance,
        replicates=replicates,seed=seed,**kwargs)
    return paths,summarize_simulations(paths)
