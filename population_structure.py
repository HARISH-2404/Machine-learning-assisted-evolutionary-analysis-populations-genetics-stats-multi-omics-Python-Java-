"""PCA-based exploratory population structure from mean-imputed SNP dosage.

This is PCA, not an admixture/STRUCTURE Bayesian ancestry estimator. Marker
columns are centered and scaled by allele-frequency standard deviation;
near-monomorphic markers are removed. Missing calls are mean-imputed.
"""
import numpy as np
import pandas as pd
from ._core import validate_dosage,save

def population_structure(genotypes,id_col="genotype_id",n_components=5,
                         population_col="population_id",output=None):
    markers=validate_dosage(genotypes,id_col,population_col if population_col in genotypes else None)
    if n_components<1: raise ValueError("n_components must be >=1.")
    x=genotypes[markers].to_numpy(float)
    means=np.nanmean(x,axis=0)
    valid=np.isfinite(means)
    x=x[:,valid]; means=means[valid]; kept=np.asarray(markers)[valid]
    x=np.where(np.isfinite(x),x,means)
    p=means/2
    sd=np.sqrt(2*p*(1-p))
    polymorphic=sd>1e-8
    x=x[:,polymorphic]; sd=sd[polymorphic]; kept=kept[polymorphic]
    if x.shape[1]==0: raise ValueError("No polymorphic markers available for PCA.")
    z=(x-means[valid][polymorphic])/sd
    u,s,vt=np.linalg.svd(z,full_matrices=False)
    k=min(n_components,u.shape[1])
    scores=u[:,:k]*s[:k]
    eig=s**2/max(1,len(genotypes)-1)
    ratio=eig/eig.sum()
    result=pd.DataFrame(scores,columns=[f"PC{i+1}" for i in range(k)])
    result.insert(0,id_col,genotypes[id_col].to_numpy())
    if population_col in genotypes:
        result.insert(1,population_col,genotypes[population_col].to_numpy())
    loadings=pd.DataFrame(vt[:k].T,columns=[f"PC{i+1}" for i in range(k)])
    loadings.insert(0,"marker_id",kept)
    variance=pd.DataFrame({"PC":[f"PC{i+1}" for i in range(k)],
                           "eigenvalue":eig[:k],"variance_explained":ratio[:k]})
    if output:
        save(result,str(output).replace(".csv","_scores.csv"))
        save(loadings,str(output).replace(".csv","_loadings.csv"))
        save(variance,str(output).replace(".csv","_variance.csv"))
    return result,loadings,variance
