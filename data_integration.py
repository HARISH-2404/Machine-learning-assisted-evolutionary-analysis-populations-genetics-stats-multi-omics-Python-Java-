"""Join filtered genotype, trial means and metadata using explicit ID checks.

Only the supplied phenotype trait is treated as a target. No simulation truth,
causal-marker flags or leakage-prone metadata are included as ML predictors.
"""
from __future__ import annotations
import pandas as pd
from ._utils import require_columns,require_unique,save_table

class DataIntegrator:
    @staticmethod
    def join_genotype_phenotype(genotype,phenotype_means,trait,
                                location=None,output_dir=None):
        require_unique(genotype,"genotype_id")
        require_columns(phenotype_means,["genotype_id","location_id",trait])
        ph=phenotype_means.copy()
        if location is not None:
            ph=ph.loc[ph.location_id==location].copy()
        if ph.empty: raise ValueError("No phenotype records match the requested location.")
        if ph.duplicated(["genotype_id","location_id"]).any():
            raise ValueError("Phenotype means contain duplicate genotype-location pairs.")
        geno_cols=[c for c in genotype if c=="genotype_id" or c.startswith("SNP_")]
        if len(geno_cols)<2:
            raise ValueError("Expected genotype_id and SNP_ marker columns.")
        joined=ph[["genotype_id","location_id",trait]].merge(
            genotype[geno_cols],on="genotype_id",how="inner",validate="many_to_one")
        if joined.empty: raise ValueError("No matching genotype IDs found.")
        # Multiple locations per genotype must be grouped in downstream CV to avoid leakage.
        if output_dir:
            save_table(joined,f"{output_dir}/integrated_{trait}.csv")
        return joined
