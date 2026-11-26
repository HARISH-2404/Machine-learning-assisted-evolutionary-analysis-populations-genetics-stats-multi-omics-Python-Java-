"""Haplotype frequencies from phased 0/1 allele columns.

Input must be phased haplotypes: rows are chromosomes/haplotypes and marker
columns contain 0 or 1. Unphased diploid dosage (0/1/2) cannot uniquely
identify haplotypes; this function intentionally does not guess phase.
"""
import pandas as pd
from ._core import save

def haplotype_frequency(phased_haplotypes,marker_cols=None,output=None):
    cols=list(marker_cols or [c for c in phased_haplotypes.columns if c not in {"haplotype_id","population_id"}])
    if not cols: raise ValueError("No phased marker columns supplied.")
    x=phased_haplotypes[cols]
    if x.isna().any().any() or not x.isin([0,1]).all().all():
        raise ValueError("Phased haplotype markers must contain only 0/1 and no missing values.")
    data=phased_haplotypes.copy()
    data["haplotype"]=x.astype(int).astype(str).agg("".join,axis=1)
    group_cols=["haplotype"]+(["population_id"] if "population_id" in data else [])
    counts=data.groupby(group_cols).size().rename("count").reset_index()
    denom=counts.groupby("population_id")["count"].transform("sum") if "population_id" in counts else counts["count"].sum()
    counts["frequency"]=counts["count"]/denom
    if output: save(counts,output)
    return counts
