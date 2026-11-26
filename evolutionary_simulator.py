"""Forward-time Wright-Fisher simulation with selection, dominance, drift,
mutation and migration.

For allele A frequency p and fitnesses:
aa=1, Aa=1+h*s, AA=1+s.
After viability selection, migration and mutation, the next generation is
sampled as Binomial(2N, p')/(2N), representing genetic drift in a diploid
Wright-Fisher population.

This is a single-locus/multi-locus independent-frequency approximation.
It does not simulate recombination, linkage disequilibrium, haplotypes,
polyploid inheritance, or epistasis.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._common import rng_from, write_csv, get_logger

class EvolutionarySimulator:
    def __init__(self, generations=50, population_size=500,
                 selection_coefficient=0.05, dominance_coefficient=0.5,
                 mutation_rate=1e-6, migration_rate=0.01,
                 initial_frequency=0.20, seed=20261127):
        if generations < 1 or population_size < 2:
            raise ValueError("generations >= 1 and population_size >= 2 are required.")
        if not 0 <= dominance_coefficient <= 1:
            raise ValueError("dominance_coefficient h must be in [0,1].")
        if not 0 <= mutation_rate <= 1 or not 0 <= migration_rate <= 1:
            raise ValueError("Mutation and migration rates must be in [0,1].")
        if not 0 <= initial_frequency <= 1:
            raise ValueError("initial_frequency must be in [0,1].")
        if selection_coefficient <= -1:
            raise ValueError("selection_coefficient must be > -1 for this fitness parameterization.")
        self.generations, self.N = generations, population_size
        self.s, self.h = selection_coefficient, dominance_coefficient
        self.mu, self.m, self.p0 = mutation_rate, migration_rate, initial_frequency
        self.rng, self.log = rng_from(seed), get_logger(self.__class__.__name__)

    @staticmethod
    def selection_step(p, s, h):
        """Return allele frequency after viability selection."""
        w_aa, w_Aa, w_AA = 1.0, 1.0+h*s, 1.0+s
        mean_w = (1-p)**2*w_aa + 2*p*(1-p)*w_Aa + p**2*w_AA
        return (p**2*w_AA + p*(1-p)*w_Aa) / np.maximum(mean_w,1e-12)

    def simulate(self, initial_frequencies=None, population_ids=None,
                 marker_ids=None, output_dir=None):
        p = np.asarray(self.p0 if initial_frequencies is None else initial_frequencies,dtype=float)
        if p.ndim == 1:
            p = p[None,:]
        if p.ndim != 2 or not np.isfinite(p).all() or ((p<0)|(p>1)).any():
            raise ValueError("initial_frequencies must be a finite 1D/2D array within [0,1].")
        n_pop,n_markers = p.shape
        pops = population_ids or [f"POP_{i+1}" for i in range(n_pop)]
        markers = marker_ids or [f"SNP_{i+1:03d}" for i in range(n_markers)]
        if len(pops)!=n_pop or len(markers)!=n_markers:
            raise ValueError("population_ids/marker_ids lengths must match input dimensions.")
        records=[]
        for generation in range(self.generations+1):
            for pi,pop in enumerate(pops):
                for mi,marker in enumerate(markers):
                    records.append({"population_id":pop,"generation":generation,
                                    "marker_id":marker,"allele_frequency":float(p[pi,mi]),
                                    "selection_coefficient":self.s,"dominance_coefficient":self.h,
                                    "population_size":self.N})
            if generation==self.generations:
                break
            selected=self.selection_step(p,self.s,self.h)
            # Symmetric island migration towards mean frequency across populations.
            migrated=(1-self.m)*selected+self.m*selected.mean(axis=0,keepdims=True)
            # Forward and reverse mutation at equal rates.
            mutated=migrated*(1-self.mu)+(1-migrated)*self.mu
            p=self.rng.binomial(2*self.N,np.clip(mutated,0,1))/(2*self.N)
        result=pd.DataFrame(records)
        if output_dir:
            write_csv(result,f"{output_dir}/allele_frequency_trajectories.csv")
            summary=(result.groupby(["population_id","marker_id"],as_index=False)
                     .agg(initial_frequency=("allele_frequency","first"),
                          final_frequency=("allele_frequency","last"),
                          min_frequency=("allele_frequency","min"),
                          max_frequency=("allele_frequency","max")))
            write_csv(summary,f"{output_dir}/evolutionary_summary.csv")
        self.log.info("Simulated %d generations for %d populations and %d markers.",self.generations,n_pop,n_markers)
        return result
