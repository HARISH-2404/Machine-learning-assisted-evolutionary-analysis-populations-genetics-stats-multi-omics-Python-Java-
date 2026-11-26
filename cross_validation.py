"""K-fold and group-aware cross-validation index generators."""
import numpy as np

def kfold_indices(n_samples,n_splits=5,shuffle=True,seed=42):
    if n_splits<2 or n_splits>n_samples: raise ValueError("Require 2 <= n_splits <= n_samples")
    idx=np.arange(n_samples)
    if shuffle: np.random.default_rng(seed).shuffle(idx)
    folds=np.array_split(idx,n_splits)
    for i in range(n_splits):
        test=folds[i]; train=np.concatenate([folds[j] for j in range(n_splits) if j!=i])
        yield train,test

def grouped_kfold_indices(groups,n_splits=5,seed=42):
    groups=np.asarray(groups); unique=np.unique(groups)
    if n_splits<2 or n_splits>len(unique): raise ValueError("n_splits must fit distinct groups")
    rng=np.random.default_rng(seed); rng.shuffle(unique); bins=np.array_split(unique,n_splits)
    for b in bins:
        test=np.where(np.isin(groups,b))[0]; train=np.where(~np.isin(groups,b))[0]
        yield train,test
