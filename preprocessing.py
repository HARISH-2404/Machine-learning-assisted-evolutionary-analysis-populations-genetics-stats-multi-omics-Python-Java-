"""Leakage-aware numeric preprocessing and deterministic data splitting."""
import numpy as np

def standardize(X, mean=None, scale=None):
    X=np.asarray(X,dtype=float)
    if X.ndim!=2 or not np.isfinite(X).all(): raise ValueError("X must be a finite 2D matrix")
    mu=X.mean(axis=0) if mean is None else np.asarray(mean,dtype=float)
    sd=X.std(axis=0) if scale is None else np.asarray(scale,dtype=float)
    if mu.shape!=(X.shape[1],) or sd.shape!=(X.shape[1],): raise ValueError("mean/scale dimensions mismatch")
    sd=sd.copy(); sd[sd==0]=1.0
    return (X-mu)/sd,mu,sd

def train_test_split_indices(n_samples, test_fraction=0.2, seed=42, groups=None):
    if n_samples<2 or not 0<test_fraction<1: raise ValueError("Need >=2 samples and test_fraction in (0,1)")
    rng=np.random.default_rng(seed)
    if groups is None:
        idx=rng.permutation(n_samples); ntest=max(1,min(n_samples-1,round(n_samples*test_fraction)))
        return idx[ntest:],idx[:ntest]
    groups=np.asarray(groups)
    if len(groups)!=n_samples: raise ValueError("groups length mismatch")
    unique=np.unique(groups); rng.shuffle(unique)
    ng=max(1,min(len(unique)-1,round(len(unique)*test_fraction)))
    test_groups=set(unique[:ng]); test=np.array([i for i,g in enumerate(groups) if g in test_groups])
    train=np.array([i for i,g in enumerate(groups) if g not in test_groups])
    if len(train)==0 or len(test)==0: raise ValueError("Need at least two distinct groups")
    return train,test
