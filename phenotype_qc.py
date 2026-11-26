"""Quality-control multi-environment phenotype tables and summarize trial data.

Supports a wide plot table with genotype_id, location_id, replication and
trait columns. Missing values are preserved by default. Outliers are flagged
using within-location robust median absolute deviation (MAD); no automatic
removal, so field investigators can review suspected biological extremes.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._utils import require_columns,save_table

class PhenotypeQC:
    def __init__(self,expected_locations=5,expected_reps=3,mad_threshold=4.0):
        if expected_locations<1 or expected_reps<1 or mad_threshold<=0:
            raise ValueError("Expected design sizes and MAD threshold must be positive.")
        self.expected_locations=expected_locations
        self.expected_reps=expected_reps
        self.mad_threshold=mad_threshold

    def inspect(self,frame,traits,output_dir=None):
        require_columns(frame,["genotype_id","location_id","replication"]+list(traits))
        if frame[["genotype_id","location_id","replication"]].isna().any().any():
            raise ValueError("Trial design identifiers cannot be missing.")
        duplicate=frame.duplicated(["genotype_id","location_id","replication"])
        if duplicate.any():
            raise ValueError(f"Duplicate genotype-location-replication plots: {int(duplicate.sum())}")
        locations=frame.location_id.nunique()
        design=(frame.groupby(["genotype_id","location_id"]).replication
                .agg(["nunique","size"]).reset_index())
        incomplete=design[(design["nunique"]!=self.expected_reps)|
                          (design["size"]!=self.expected_reps)]
        if locations!=self.expected_locations:
            raise ValueError(f"Expected {self.expected_locations} locations; found {locations}.")
        # Detect missing genotype-location combinations, not just incomplete observed groups.
        expected_pairs=frame.genotype_id.nunique()*self.expected_locations
        if len(design)!=expected_pairs:
            raise ValueError("Some genotype-location combinations are absent.")
        if not incomplete.empty:
            raise ValueError(f"{len(incomplete)} genotype-location combinations lack required reps.")
        flagged=frame.copy()
        summaries=[]
        for trait in traits:
            values=pd.to_numeric(frame[trait],errors="coerce")
            flagged[trait]=values
            outlier=np.zeros(len(frame),dtype=bool)
            for _,idx in frame.groupby("location_id").groups.items():
                v=values.loc[idx]
                med=v.median()
                mad=(v-med).abs().median()
                if pd.notna(mad) and mad>0:
                    outlier[frame.index.get_indexer(idx)]=(0.67448975*(v-med).abs()/mad>self.mad_threshold).fillna(False)
            flagged[f"{trait}_outlier_flag"]=outlier
            summaries.append({"trait":trait,"observations":int(values.notna().sum()),
                              "missing":int(values.isna().sum()),
                              "outliers_flagged":int(outlier.sum()),
                              "mean":values.mean(),"std":values.std()})
        summary=pd.DataFrame(summaries)
        means=(flagged.groupby(["genotype_id","location_id"],as_index=False)
               [list(traits)].mean())
        if output_dir:
            save_table(flagged,f"{output_dir}/phenotypes_with_qc_flags.csv")
            save_table(summary,f"{output_dir}/phenotype_trait_summary.csv")
            save_table(means,f"{output_dir}/genotype_location_means.csv")
        return flagged,summary,means
