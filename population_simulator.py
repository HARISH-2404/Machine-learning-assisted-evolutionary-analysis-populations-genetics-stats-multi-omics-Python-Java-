"""Create population sample metadata and allele-frequency summaries."""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._common import rng_from, write_csv, validate_diploid_dosage, get_logger

class PopulationSimulator:
    def __init__(self, n_populations=5, individuals_per_population=50,
                 generations=10, population_names=None, seed=20261127):
        if min(n_populations, individuals_per_population, generations) < 1:
            raise ValueError("All population dimensions must be >= 1.")
        self.n_populations, self.n_individuals = n_populations, individuals_per_population
        self.generations = generations
        self.names = population_names or [f"POP_{i+1}" for i in range(n_populations)]
        if len(self.names) != n_populations or len(set(self.names)) != len(self.names):
            raise ValueError("Population names must be unique and match n_populations.")
        self.rng, self.log = rng_from(seed), get_logger(self.__class__.__name__)

    def simulate(self, n_markers=300, genotype=None, output_dir=None):
        """Return individual metadata and marker-frequency table.

        If genotype is supplied it must contain population_id in its metadata
        columns; frequencies are calculated from the corresponding dosage rows.
        """
        if n_markers < 1:
            raise ValueError("n_markers must be positive.")
        sample_rows = []
        for pop in self.names:
            for i in range(1,self.n_individuals+1):
                sample_rows.append({"sample_id":f"{pop}_I{i:04d}",
                                    "population_id":pop,"generation":0,
                                    "individual_index":i,"data_origin":"synthetic"})
        samples = pd.DataFrame(sample_rows)
        freq_rows = []
        if genotype is not None:
            if "population_id" not in genotype:
                raise ValueError("Genotype data must include population_id for population frequency summaries.")
            # Validate only genotype_id + marker columns; population_id is metadata.
            marker_frame = genotype.drop(columns=["population_id"])
            marker_cols = validate_diploid_dosage(marker_frame)
            for pop in self.names:
                subset = genotype.loc[genotype.population_id == pop, marker_cols].astype(float)
                if subset.empty:
                    continue
                for marker in marker_cols:
                    called = subset[marker].dropna()
                    p = called.sum()/(2*len(called)) if len(called) else np.nan
                    freq_rows.append({"population_id":pop,"generation":0,
                                      "marker_id":marker,"allele_frequency":p,
                                      "sample_size":len(called)})
        else:
            p0 = self.rng.uniform(0.05,0.5,(len(self.names),n_markers))
            for pi,pop in enumerate(self.names):
                for mi in range(n_markers):
                    freq_rows.append({"population_id":pop,"generation":0,
                                      "marker_id":f"SNP_{mi+1:03d}",
                                      "allele_frequency":p0[pi,mi],
                                      "sample_size":self.n_individuals})
        frequencies = pd.DataFrame(freq_rows)
        if output_dir:
            write_csv(samples,f"{output_dir}/population_membership.csv")
            write_csv(frequencies,f"{output_dir}/population_allele_frequencies.csv")
        self.log.info("Generated metadata for %d individuals across %d populations.",len(samples),len(self.names))
        return samples, frequencies
