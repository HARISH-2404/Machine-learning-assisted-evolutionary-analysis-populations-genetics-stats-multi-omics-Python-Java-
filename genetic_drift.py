"""Neutral Wright-Fisher genetic drift for diploid populations."""
import numpy as np
from ._core import prob,ne_int

def simulate_drift(initial_frequency,generations,effective_population_size,replicates=1,seed=None):
    """Return [replicate,generation] frequencies; each step samples 2*Ne alleles."""
    p=prob(initial_frequency,'initial_frequency'); ne=ne_int(effective_population_size)
    if int(generations)!=generations or generations<0 or int(replicates)!=replicates or replicates<1: raise ValueError('invalid generations or replicates')
    rng=np.random.default_rng(seed); out=np.empty((int(replicates),int(generations)+1)); out[:,0]=p
    for r in range(int(replicates)):
        for t in range(1,int(generations)+1): out[r,t]=rng.binomial(2*ne,out[r,t-1])/(2*ne)
    return out
