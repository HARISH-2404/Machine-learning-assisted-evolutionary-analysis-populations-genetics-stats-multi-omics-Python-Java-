"""Simulate diploid SNP dosage, marker map, QC and population structure.

The model uses Hardy-Weinberg genotype sampling at each marker:
P(AA)=(1-p)^2, P(Aa)=2p(1-p), P(aa)=p^2.
This is a teaching/demo model; it does not simulate linkage disequilibrium
unless an explicit shared-haplotype extension is added.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._common import rng_from, write_csv, validate_diploid_dosage, summarize_dosage, get_logger

class GenotypeSimulator:
    def __init__(self, n_genotypes=250, n_markers=300, n_chromosomes=15,
                 maf_min=0.05, maf_max=0.50, missing_rate=0.01,
                 n_populations=5, population_divergence=0.12, seed=20261127):
        if min(n_genotypes, n_markers, n_chromosomes, n_populations) < 1:
            raise ValueError("All dimensions must be positive integers.")
        if not 0 <= maf_min <= maf_max <= 0.5:
            raise ValueError("Require 0 <= maf_min <= maf_max <= 0.5.")
        if not 0 <= missing_rate < 1:
            raise ValueError("missing_rate must be in [0, 1).")
        if not 0 <= population_divergence <= 0.5:
            raise ValueError("population_divergence must be in [0, 0.5].")
        self.n_genotypes, self.n_markers = n_genotypes, n_markers
        self.n_chromosomes, self.maf_min, self.maf_max = n_chromosomes, maf_min, maf_max
        self.missing_rate, self.n_populations = missing_rate, n_populations
        self.population_divergence, self.rng = population_divergence, rng_from(seed)
        self.log = get_logger(self.__class__.__name__)

    def simulate(self, output_dir=None, prefix="simulated"):
        """Return genotype dosage, marker map and sample metadata."""
        ids = [f"G{i:03d}" for i in range(1, self.n_genotypes+1)]
        marker_ids = [f"SNP_{i:03d}" for i in range(1, self.n_markers+1)]
        pop_ids = [f"POP_{(i % self.n_populations)+1}" for i in range(self.n_genotypes)]
        base_p = self.rng.uniform(self.maf_min, self.maf_max, self.n_markers)
        pop_shift = self.rng.normal(0, self.population_divergence, (self.n_populations, self.n_markers))
        pop_freq = np.clip(base_p[None, :] + pop_shift, 0.001, 0.999)
        dosage = np.empty((self.n_genotypes, self.n_markers), dtype=float)
        for i, pop in enumerate(pop_ids):
            p = pop_freq[int(pop.split("_")[1])-1]
            dosage[i] = self.rng.binomial(2, p)
        dosage[self.rng.random(dosage.shape) < self.missing_rate] = np.nan
        genotype = pd.DataFrame(dosage, columns=marker_ids)
        genotype.insert(0, "genotype_id", ids)
        per_chr = int(np.ceil(self.n_markers/self.n_chromosomes))
        marker_map = pd.DataFrame({
            "marker_id": marker_ids,
            "chromosome": [f"Chr{min(i//per_chr+1, self.n_chromosomes):02d}" for i in range(self.n_markers)],
            "position_bp": [(i%per_chr+1)*100_000 for i in range(self.n_markers)],
            "reference_allele": "A", "alternate_allele": "G",
            "simulated_source_maf": base_p, "marker_type": "SNP",
            "simulation_truth_causal": False
        })
        metadata = pd.DataFrame({"genotype_id": ids, "population_id": pop_ids,
                                 "ploidy_assumption": 2, "data_origin": "synthetic"})
        validate_diploid_dosage(genotype)
        if output_dir:
            write_csv(genotype, f"{output_dir}/{prefix}_genotype_dosage.csv")
            write_csv(marker_map, f"{output_dir}/{prefix}_marker_map.csv")
            write_csv(metadata, f"{output_dir}/{prefix}_sample_metadata.csv")
            write_csv(summarize_dosage(genotype), f"{output_dir}/{prefix}_marker_qc.csv")
        self.log.info("Generated %d genotypes x %d markers.", self.n_genotypes, self.n_markers)
        return genotype, marker_map, metadata

    @staticmethod
    def validate(genotype: pd.DataFrame) -> list[str]:
        return validate_diploid_dosage(genotype)
