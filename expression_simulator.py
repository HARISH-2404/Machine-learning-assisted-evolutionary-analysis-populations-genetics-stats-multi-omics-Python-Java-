"""Generate synthetic RNA-seq-like raw counts using a negative-binomial model.

Counts are simulated as NB(mean=mu, dispersion=alpha), where
Var(Y)=mu + alpha*mu^2. These are toy counts for pipeline testing only.
They are not biological replicates or real RNA-seq evidence.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from ._common import rng_from, write_csv, get_logger

class ExpressionSimulator:
    def __init__(self, n_genes=100, mean_library_size=2_000_000,
                 dispersion=0.10, de_fraction=0.10, log2_fold_change=1.0,
                 seed=20261127):
        if n_genes < 1 or mean_library_size < 1 or dispersion <= 0:
            raise ValueError("n_genes and mean_library_size must be positive; dispersion > 0.")
        if not 0 <= de_fraction <= 1:
            raise ValueError("de_fraction must be in [0,1].")
        self.n_genes, self.mean_library_size = n_genes, mean_library_size
        self.dispersion, self.de_fraction = dispersion, de_fraction
        self.log2_fold_change, self.rng = log2_fold_change, rng_from(seed)
        self.log = get_logger(self.__class__.__name__)

    def simulate(self, sample_ids, groups=None, output_dir=None):
        sample_ids = list(map(str, sample_ids))
        n = len(sample_ids)
        if not n or len(set(sample_ids)) != n:
            raise ValueError("sample_ids must be non-empty and unique.")
        groups = np.asarray(groups if groups is not None else ["control"]*n, dtype=str)
        if len(groups) != n:
            raise ValueError("groups must contain one label per sample.")
        levels = list(dict.fromkeys(groups.tolist()))
        genes = [f"GENE_{i:04d}" for i in range(1,self.n_genes+1)]
        base_mean = self.rng.lognormal(mean=3.0, sigma=1.0, size=self.n_genes)
        de_n = int(round(self.n_genes*self.de_fraction))
        de_idx = set(self.rng.choice(self.n_genes, de_n, replace=False).tolist()) if de_n else set()
        lib_factors = self.rng.lognormal(0, 0.15, n)
        counts = np.zeros((n,self.n_genes), dtype=np.int64)
        truth = []
        for j, gene in enumerate(genes):
            truth.append({"gene_id":gene, "simulated_differential":j in de_idx,
                          "simulated_log2_fold_change":self.log2_fold_change if j in de_idx else 0.0})
            for i in range(n):
                fold = 2**self.log2_fold_change if len(levels)>1 and groups[i] != levels[0] and j in de_idx else 1.0
                mu = max(0.01, base_mean[j]*lib_factors[i]*fold)
                size = 1/self.dispersion
                prob = size/(size+mu)
                counts[i,j] = self.rng.negative_binomial(size, prob)
        result = pd.DataFrame(counts, columns=genes)
        result.insert(0,"group",groups)
        result.insert(0,"sample_id",sample_ids)
        truth_df = pd.DataFrame(truth)
        if output_dir:
            write_csv(result, f"{output_dir}/expression_counts.csv")
            write_csv(truth_df, f"{output_dir}/expression_simulation_truth.csv")
        self.log.info("Generated %d samples x %d genes.", n, self.n_genes)
        return result, truth_df
