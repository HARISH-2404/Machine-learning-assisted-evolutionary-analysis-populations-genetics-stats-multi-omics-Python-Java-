"""Evidence-aware tracking of candidate resistance-associated alleles."""
import pandas as pd

def track_resistance_candidate(trajectory, allele_label, validated_effect=False,
                               threshold=0.8, group_cols=("population","location")):
    """Label allele frequency status; never infer resistance from frequency alone."""
    x=trajectory[trajectory["allele"].astype(str)==str(allele_label)].copy()
    if x.empty: raise ValueError("Requested allele not found")
    x["candidate_allele_label"]=("validated_resistance_associated" if validated_effect
                                  else "unvalidated_candidate_allele")
    x["frequency_status"]=x.allele_frequency.map(
        lambda p:"high_frequency" if p>=threshold else "below_threshold")
    x["interpretation"]="frequency trajectory only; resistance phenotype requires independent validation"
    return x.sort_values([*group_cols,"generation"])

def compare_resistance_and_susceptibility(trajectory, resistance_allele, susceptible_allele,
                                         group_cols=("population","location","generation")):
    x=trajectory[trajectory.allele.astype(str).isin([str(resistance_allele),str(susceptible_allele)])]
    wide=x.pivot_table(index=list(group_cols),columns="allele",values="allele_frequency",aggfunc="mean").reset_index()
    if resistance_allele not in wide or susceptible_allele not in wide:
        raise ValueError("Both allele classes must be present")
    wide["resistance_minus_susceptible_frequency"]=wide[resistance_allele]-wide[susceptible_allele]
    return wide
