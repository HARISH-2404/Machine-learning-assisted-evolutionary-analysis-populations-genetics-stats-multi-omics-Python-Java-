"""Longitudinal locus-to-population-to-trait integration.

This module summarizes allele-frequency change at a genomic locus across
populations, locations, sampling times and generations, and joins optional
quantitative-genetic/phenotypic observations. It does not infer causality.
"""
import pandas as pd
import numpy as np

REQUIRED={"chromosome","position","gene_id","allele","population","location",
          "timepoint","generation","allele_count","total_alleles"}

def build_evolutionary_trajectory(allele_counts, trait_data=None, gene_annotations=None):
    """Build locus x population x place x time x generation trajectory table.

    allele_counts: one row per allele/locus/group, with allele_count and
    total_alleles. Frequencies are calculated as allele_count/total_alleles.
    Optional trait_data should share population, location, timepoint, generation.
    """
    missing=REQUIRED-set(allele_counts.columns)
    if missing: raise ValueError(f"Missing required columns: {sorted(missing)}")
    x=allele_counts.copy()
    if (x["total_alleles"]<=0).any() or (x["allele_count"]<0).any() or (x["allele_count"]>x["total_alleles"]).any():
        raise ValueError("Allele counts must satisfy 0 <= allele_count <= total_alleles and total_alleles > 0")
    x["allele_frequency"]=x["allele_count"]/x["total_alleles"]
    keys=["chromosome","position","gene_id","allele","population","location","timepoint","generation"]
    x=x.sort_values(["chromosome","position","population","location","generation","timepoint"])
    grp=["chromosome","position","gene_id","allele","population","location"]
    x["previous_frequency"]=x.groupby(grp,dropna=False)["allele_frequency"].shift(1)
    x["frequency_change"]=x["allele_frequency"]-x["previous_frequency"]
    x["frequency_change_per_generation"]=x["frequency_change"] / (
        x["generation"]-x.groupby(grp,dropna=False)["generation"].shift(1)).replace(0,np.nan)
    if gene_annotations is not None:
        x=x.merge(gene_annotations,on=["gene_id"],how="left",suffixes=("","_annotation"))
    if trait_data is not None:
        join_keys=["population","location","timepoint","generation"]
        absent=set(join_keys)-set(trait_data.columns)
        if absent: raise ValueError(f"trait_data missing join keys: {sorted(absent)}")
        x=x.merge(trait_data,on=join_keys,how="left",validate="many_to_many")
    return x

def classify_allele_state(trajectory, resistance_allele, susceptible_allele=None,
                          frequency_threshold=0.8):
    """Label frequency pattern, not proven biological resistance causality.

    Requires an explicit allele label known from independent evidence. A rise
    in frequency alone is not sufficient to call an allele resistance-conferring.
    """
    x=trajectory.copy()
    if "allele" not in x or "allele_frequency" not in x:
        raise ValueError("trajectory requires allele and allele_frequency")
    x["state_label"]="other_or_unclassified"
    x.loc[x["allele"].astype(str)==str(resistance_allele),"state_label"]="candidate_resistance_allele"
    if susceptible_allele is not None:
        x.loc[x["allele"].astype(str)==str(susceptible_allele),"state_label"]="candidate_susceptible_allele"
    x["high_frequency_flag"]=x["allele_frequency"]>=frequency_threshold
    return x

def summarize_locus_population(trajectory):
    """Summarize latest observed frequency and net change per locus/population/place."""
    keys=["chromosome","position","gene_id","allele","population","location"]
    x=trajectory.sort_values(["generation","timepoint"])
    first=x.groupby(keys,dropna=False).first().reset_index()
    last=x.groupby(keys,dropna=False).last().reset_index()
    out=last[keys+["timepoint","generation","allele_frequency"]].rename(
        columns={"timepoint":"latest_timepoint","generation":"latest_generation",
                 "allele_frequency":"latest_frequency"})
    out["initial_frequency"]=first["allele_frequency"].to_numpy()
    out["net_frequency_change"]=out["latest_frequency"]-out["initial_frequency"]
    return out
