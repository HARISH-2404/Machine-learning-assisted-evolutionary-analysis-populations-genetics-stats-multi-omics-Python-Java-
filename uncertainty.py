"""Monte Carlo summaries and probability-of-threshold calculations."""
import numpy as np
import pandas as pd

def probability_above(paths, threshold, at_generation=None):
    x=paths if at_generation is None else paths[paths.generation==at_generation]
    if x.empty: raise ValueError("No simulation rows for requested generation")
    return float((x.allele_frequency>=threshold).mean())

def fixation_loss_probability(paths, at_generation=None, tolerance=0.0):
    x=paths if at_generation is None else paths[paths.generation==at_generation]
    if x.empty: raise ValueError("No simulation rows for requested generation")
    return {"fixation_probability":float((x.allele_frequency>=1-tolerance).mean()),
            "loss_probability":float((x.allele_frequency<=tolerance).mean())}

def bootstrap_mean_interval(values, confidence=0.95, n_boot=2000, seed=1):
    x=np.asarray(values,dtype=float); x=x[np.isfinite(x)]
    if len(x)<2: raise ValueError("At least two finite values required")
    rng=np.random.default_rng(seed)
    means=np.mean(rng.choice(x,(n_boot,len(x)),replace=True),axis=1)
    a=(1-confidence)/2
    return {"mean":float(x.mean()),"lower":float(np.quantile(means,a)),
            "upper":float(np.quantile(means,1-a))}
