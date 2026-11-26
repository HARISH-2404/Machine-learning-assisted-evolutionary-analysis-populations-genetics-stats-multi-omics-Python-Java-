"""Wright-Fisher forward simulation with selection, drift, mutation and migration."""
import numpy as np
import pandas as pd
from .fitness_models import genotype_fitness

def deterministic_diploid_selection(p, fitnesses):
    """Allele frequency after viability selection under random mating."""
    if not 0<=p<=1: raise ValueError("p must be in [0,1]")
    wAA,wAa,waa=fitnesses; q=1-p
    mean=p*p*wAA+2*p*q*wAa+q*q*waa
    if mean<=0: raise ValueError("Mean fitness must be positive")
    return (p*p*wAA+p*q*wAa)/mean

def simulate_wright_fisher(p0, generations, population_size, selection_coefficient=0.0,
                           dominance=0.5, fitness_model="additive", mutation_forward=0.0,
                           mutation_reverse=0.0, migration_rate=0.0, migrant_frequency=None,
                           replicates=100, seed=1, fitnesses=None, environment_effects=None):
    """Simulate diploid biallelic allele-frequency paths.

    Drift is binomial sampling of 2N gene copies after deterministic selection,
    mutation and migration. `environment_effects` optionally maps generation to
    an additive adjustment to s. The model assumes a single well-defined locus.
    """
    if not 0<=p0<=1 or generations<0 or population_size<1 or replicates<1:
        raise ValueError("Invalid p0, generations, population_size or replicates")
    if not (0<=mutation_forward<=1 and 0<=mutation_reverse<=1 and 0<=migration_rate<=1):
        raise ValueError("Mutation and migration rates must be in [0,1]")
    if migration_rate and migrant_frequency is None: raise ValueError("migrant_frequency required when migration_rate > 0")
    if migrant_frequency is not None and not 0<=migrant_frequency<=1: raise ValueError("migrant_frequency must be in [0,1]")
    rng=np.random.default_rng(seed); rows=[]
    for rep in range(replicates):
        p=float(p0)
        rows.append((rep,0,p))
        for gen in range(1,generations+1):
            s=selection_coefficient+(environment_effects or {}).get(gen,0.0)
            fit=fitnesses if fitnesses is not None else genotype_fitness(s,dominance,fitness_model)
            p=deterministic_diploid_selection(p,fit)
            p=p*(1-mutation_forward)+(1-p)*mutation_reverse
            if migration_rate:
                p=(1-migration_rate)*p+migration_rate*migrant_frequency
            copies=rng.binomial(2*population_size,p)
            p=copies/(2*population_size)
            rows.append((rep,gen,p))
    return pd.DataFrame(rows,columns=["replicate","generation","allele_frequency"])

def summarize_simulations(paths, confidence=0.95):
    if not 0<confidence<1: raise ValueError("confidence must be between 0 and 1")
    alpha=(1-confidence)/2
    g=paths.groupby("generation")["allele_frequency"]
    out=g.agg(mean="mean",sd="std",median="median").reset_index()
    out["lower_quantile"]=g.quantile(alpha).to_numpy()
    out["upper_quantile"]=g.quantile(1-alpha).to_numpy()
    return out
