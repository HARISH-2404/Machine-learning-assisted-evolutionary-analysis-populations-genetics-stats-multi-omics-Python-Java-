"""Shared numerical validation and diploid biallelic equations."""
import numpy as np

def prob(x, name='probability'):
    x=float(x)
    if not 0 <= x <= 1: raise ValueError(f'{name} must be in [0,1]')
    return x

def ne_int(x):
    if int(x)!=x or int(x)<1: raise ValueError('effective population size must be a positive integer')
    return int(x)

def genotype_frequencies(p):
    p=prob(p,'allele frequency'); q=1-p
    return np.array([p*p,2*p*q,q*q],dtype=float)  # AA,Aa,aa

def selected_frequency(p, fitness):
    f=np.asarray(fitness,dtype=float)
    if f.shape!=(3,) or np.any(f<0) or not np.all(np.isfinite(f)): raise ValueError('fitness must be [wAA,wAa,waa], finite and nonnegative')
    g=genotype_frequencies(p); mean=float(g@f)
    if mean==0: raise ValueError('mean fitness is zero')
    return float((g[0]*f[0]+0.5*g[1]*f[1])/mean),mean
