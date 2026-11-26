"""Breeding response and genetic gain projections across cycles."""
import pandas as pd

def breeder_equation_response(heritability, selection_differential):
    """R = h2 * S; h2 is narrow-sense heritability on matching scale."""
    if not 0<=heritability<=1: raise ValueError("heritability must be in [0,1]")
    return float(heritability*selection_differential)

def project_breeding_response(initial_mean, generations, heritability,
                              selection_differential, genetic_gain_per_cycle=None):
    if generations<0: raise ValueError("generations must be nonnegative")
    gain=breeder_equation_response(heritability,selection_differential) if genetic_gain_per_cycle is None else float(genetic_gain_per_cycle)
    return pd.DataFrame({"generation":range(generations+1),
                         "predicted_trait_mean":[initial_mean+i*gain for i in range(generations+1)],
                         "genetic_gain_from_baseline":[i*gain for i in range(generations+1)]})

def selection_response_from_observations(parent_mean, offspring_mean, generations=1):
    if generations<=0: raise ValueError("generations must be positive")
    response=float(offspring_mean-parent_mean)
    return {"realized_response":response,"response_per_generation":response/generations}
