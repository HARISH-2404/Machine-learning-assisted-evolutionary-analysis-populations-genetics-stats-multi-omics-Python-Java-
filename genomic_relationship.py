"""Genomic relationship matrices using VanRaden-style marker scaling."""
import numpy as np

def vanraden_relationship(genotypes, allele_frequencies=None, method=1):
    """Compute additive genomic relationship matrix from diploid dosage.

    Parameters
    ----------
    genotypes : ndarray (n_individuals, n_markers), values 0/1/2; NaN allowed.
    allele_frequencies : optional ALT frequencies for each marker.
    method : 1 uses Z Z' / sum(2p(1-p)); 2 uses marker-wise standardized Z.
    Missing genotypes are mean-imputed at the marker using observed dosage.
    """
    G=np.asarray(genotypes,dtype=float)
    if G.ndim!=2: raise ValueError("genotypes must be individuals x markers")
    n,m=G.shape
    if np.any((G[np.isfinite(G)]<0)|(G[np.isfinite(G)]>2)): raise ValueError("dosages must be in [0,2]")
    if allele_frequencies is None:
        p=np.nanmean(G,axis=0)/2
    else:
        p=np.asarray(allele_frequencies,dtype=float)
        if p.shape!=(m,) or np.any((p<0)|(p>1)): raise ValueError("allele_frequencies must have one value in [0,1] per marker")
    valid=np.isfinite(p)&(p>0)&(p<1)
    G=G[:,valid]; p=p[valid]
    if G.shape[1]==0: raise ValueError("no polymorphic markers available")
    means=2*p
    G=G.copy()
    for j, mean in enumerate(means):
        missing=~np.isfinite(G[:,j])
        G[missing,j]=mean
    Z=G-2*p
    if method==1:
        denom=float(np.sum(2*p*(1-p)))
        if denom<=0: raise ValueError("relationship denominator is zero")
        K=Z@Z.T/denom
    elif method==2:
        Z=Z/np.sqrt(2*p*(1-p))
        K=Z@Z.T/Z.shape[1]
    else: raise ValueError("method must be 1 or 2")
    return (K+K.T)/2

def relationship_summary(K, individual_ids=None):
    """Return a tidy long-form table of pairwise genomic relationships."""
    K=np.asarray(K,dtype=float)
    if K.ndim!=2 or K.shape[0]!=K.shape[1]: raise ValueError("K must be square")
    ids=list(individual_ids) if individual_ids is not None else [f"IND{i+1}" for i in range(len(K))]
    if len(ids)!=len(K): raise ValueError("individual_ids length mismatch")
    return [{"individual_1":ids[i],"individual_2":ids[j],"relationship":float(K[i,j])}
            for i in range(len(K)) for j in range(i,len(K))]
