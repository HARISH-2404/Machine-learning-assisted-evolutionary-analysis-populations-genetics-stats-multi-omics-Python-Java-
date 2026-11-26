"""Fitness parameterizations for diploid genotype classes."""
import numpy as np
from ._core import prob

def genotype_fitness(selection_coefficient, dominance=0.5, baseline=1.0):
    """Return [wAA,wAa,waa] = baseline*[1+s, 1+h*s, 1].

    h=0: beneficial allele completely recessive; h=1: completely dominant;
    h=0.5: additive. This parameterization does not describe overdominance.
    """
    s=float(selection_coefficient); h=prob(dominance,'dominance'); b=float(baseline)
    if s < -1 or b <= 0: raise ValueError('selection_coefficient must be >= -1 and baseline > 0')
    return b*np.array([1+s,1+h*s,1.0])

def overdominance_fitness(w_AA,w_Aa,w_aa):
    """Explicit heterozygote-advantage fitness triplet."""
    f=np.asarray([w_AA,w_Aa,w_aa],float)
    if np.any(f<0) or not np.all(np.isfinite(f)) or not (f[1]>f[0] and f[1]>f[2]): raise ValueError('require finite nonnegative fitness and wAa > both homozygotes')
    return f

def underdominance_fitness(w_AA,w_Aa,w_aa):
    """Explicit heterozygote-disadvantage fitness triplet."""
    f=np.asarray([w_AA,w_Aa,w_aa],float)
    if np.any(f<0) or not np.all(np.isfinite(f)) or not (f[1]<f[0] and f[1]<f[2]): raise ValueError('require finite nonnegative fitness and wAa < both homozygotes')
    return f
