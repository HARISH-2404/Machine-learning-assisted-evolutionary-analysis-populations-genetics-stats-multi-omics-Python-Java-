"""Validate IDs and experimental design metadata across input tables."""
from __future__ import annotations
import pandas as pd
from ._utils import require_columns,require_unique,save_table

class MetadataQC:
    def __init__(self,expected_genotypes=250,expected_locations=5,expected_reps=3):
        if min(expected_genotypes,expected_locations,expected_reps)<1:
            raise ValueError("Expected counts must be positive.")
        self.expected_genotypes=expected_genotypes
        self.expected_locations=expected_locations
        self.expected_reps=expected_reps

    def validate_design(self,frame,output_dir=None):
        keys=["genotype_id","location_id","replication"]
        require_columns(frame,keys)
        if frame[keys].isna().any().any():
            raise ValueError("Design keys must not contain missing values.")
        if frame.duplicated(keys).any():
            raise ValueError("Repeated genotype-location-replication combinations detected.")
        gids=frame.genotype_id.nunique()
        locs=frame.location_id.nunique()
        expected=self.expected_genotypes*self.expected_locations*self.expected_reps
        if gids!=self.expected_genotypes or locs!=self.expected_locations:
            raise ValueError(f"Expected {self.expected_genotypes} genotypes and {self.expected_locations} locations; got {gids} and {locs}.")
        if len(frame)!=expected:
            raise ValueError(f"Expected {expected} plots; got {len(frame)}.")
        counts=frame.groupby(["genotype_id","location_id"]).replication.nunique()
        if len(counts)!=gids*locs or not counts.eq(self.expected_reps).all():
            raise ValueError("Unbalanced design or missing genotype-location combinations.")
        summary=pd.DataFrame([{"n_genotypes":gids,"n_locations":locs,
                               "replications":self.expected_reps,"n_plots":len(frame),
                               "balanced":True}])
        if output_dir:
            save_table(summary,f"{output_dir}/design_validation_summary.csv")
        return summary

    @staticmethod
    def check_genotype_ids(genotype,phenotype):
        require_unique(genotype,"genotype_id")
        require_columns(phenotype,["genotype_id"])
        a=set(genotype.genotype_id.astype(str))
        b=set(phenotype.genotype_id.astype(str))
        return {"only_in_genotype":sorted(a-b),
                "only_in_phenotype":sorted(b-a),
                "shared":len(a&b)}
