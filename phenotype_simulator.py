"""Simulate quantitative traits in replicated multi-environment trials.

For genotype i, environment j and replicate k:
Y_ijk = mu + G_i + E_j + (GxE)_ij + Block_jk + error_ijk.
The model also allows an additive marker score and a dominance marker score.
The trait values and units are illustrative and must be replaced for real studies.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._common import rng_from, write_csv, validate_diploid_dosage, get_logger

DEFAULT_TRAITS = {
    "seed_cotton_yield_kg_plot": (2.5, 0.55),
    "boll_weight_g": (4.8, 0.45),
    "fiber_length_mm": (29.0, 1.8),
    "fiber_strength_g_tex": (29.0, 2.0),
    "micronaire": (4.2, 0.35),
    "lint_percent": (35.0, 2.0),
    "disease_score": (3.0, 0.8),
}

class PhenotypeSimulator:
    def __init__(self, n_locations=5, replications=3, traits=None,
                 additive_fraction=0.60, dominance_fraction=0.20,
                 gxe_sd=0.25, residual_sd=0.30, block_sd=0.08, seed=20261127):
        if n_locations < 1 or replications < 1:
            raise ValueError("n_locations and replications must be >= 1.")
        if additive_fraction < 0 or dominance_fraction < 0:
            raise ValueError("Genetic variance fractions cannot be negative.")
        if additive_fraction + dominance_fraction > 1:
            raise ValueError("Additive and dominance fractions cannot sum to > 1.")
        if min(gxe_sd, residual_sd, block_sd) < 0:
            raise ValueError("Standard deviations must be non-negative.")
        self.n_locations, self.replications = n_locations, replications
        self.traits = traits or DEFAULT_TRAITS
        self.additive_fraction, self.dominance_fraction = additive_fraction, dominance_fraction
        self.gxe_sd, self.residual_sd, self.block_sd = gxe_sd, residual_sd, block_sd
        self.rng, self.log = rng_from(seed), get_logger(self.__class__.__name__)

    def simulate(self, genotype, marker_map=None, output_dir=None):
        markers = validate_diploid_dosage(genotype)
        x = genotype[markers].astype(float)
        x = x.fillna(x.mean().fillna(1.0))
        centered = x.to_numpy() - x.to_numpy().mean(axis=0, keepdims=True)
        hetero = (x.to_numpy() == 1).astype(float)
        hetero -= hetero.mean(axis=0, keepdims=True)
        gids = genotype.genotype_id.astype(str).tolist()
        locs = [f"LOC_{i}" for i in range(1, self.n_locations+1)]
        rows = []
        for trait, (mu, scale) in self.traits.items():
            additive_weights = self.rng.normal(size=len(markers))
            dominance_weights = self.rng.normal(size=len(markers))
            additive_weights /= np.linalg.norm(additive_weights) or 1
            dominance_weights /= np.linalg.norm(dominance_weights) or 1
            additive = centered @ additive_weights * scale * np.sqrt(self.additive_fraction)
            dominance = hetero @ dominance_weights * scale * np.sqrt(self.dominance_fraction)
            env = dict(zip(locs, self.rng.normal(0, 0.20*scale, len(locs))))
            gxe_matrix = self.rng.normal(0, self.gxe_sd*scale, (len(gids), len(locs)))
            for gi, gid in enumerate(gids):
                for li, loc in enumerate(locs):
                    for rep in range(1, self.replications+1):
                        block = self.rng.normal(0, self.block_sd*scale)
                        residual = self.rng.normal(0, self.residual_sd*scale)
                        value = mu + additive[gi] + dominance[gi] + env[loc] + gxe_matrix[gi,li] + block + residual
                        if trait == "disease_score":
                            value = float(np.clip(value, 0.1, 9.0))
                        else:
                            value = max(0.01, float(value))
                        rows.append({"plot_id": f"{gid}_{loc}_R{rep}", "genotype_id": gid,
                                     "location_id": loc, "replication": rep,
                                     "block_id": f"{loc}_R{rep}", "trait": trait,
                                     "value": value, "unit": "illustrative"})
        long = pd.DataFrame(rows)
        wide = long.pivot(index=["plot_id","genotype_id","location_id","replication","block_id"],
                          columns="trait", values="value").reset_index()
        means = long.groupby(["genotype_id","location_id","trait"], as_index=False).value.mean()
        if output_dir:
            write_csv(long, f"{output_dir}/phenotypes_long.csv")
            write_csv(wide, f"{output_dir}/phenotypes_wide.csv")
            write_csv(means, f"{output_dir}/genotype_location_trait_means.csv")
        expected = len(gids)*len(locs)*self.replications*len(self.traits)
        if len(long) != expected:
            raise RuntimeError(f"Expected {expected} trait observations; got {len(long)}.")
        self.log.info("Generated %d plot-trait observations.", len(long))
        return long
