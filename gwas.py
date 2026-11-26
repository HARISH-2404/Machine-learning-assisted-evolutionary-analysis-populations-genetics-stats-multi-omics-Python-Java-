"""Genome-wide association scans using single-marker regression.

This is a transparent educational implementation for quantitative traits.
It supports optional covariates and principal components. For large studies,
use a validated mixed-model GWAS package to control population structure and
relatedness (e.g. GEMMA, GAPIT, PLINK, or BOLT-LMM as appropriate).
"""
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

def _design_matrix(covariates, n):
    if covariates is None: return np.ones((n,1),dtype=float)
    c=np.asarray(covariates,dtype=float)
    if c.ndim==1: c=c[:,None]
    if c.shape[0]!=n: raise ValueError("covariates rows must match samples")
    if not np.all(np.isfinite(c)): raise ValueError("covariates contain missing/non-finite values")
    return np.column_stack([np.ones(n),c])

def single_marker_gwas(genotypes, phenotype, covariates=None, min_samples=10):
    """Fit y = intercept + covariates + SNP dosage for each marker.

    Parameters
    ----------
    genotypes : array-like, shape (n_samples, n_markers)
        Diploid biallelic alternate dosage 0/1/2; NaN allowed.
    phenotype : array-like, shape (n_samples,)
        Continuous trait values; missing phenotype samples are omitted per SNP.
    covariates : optional array-like, shape (n_samples, n_covariates)
        Include e.g. sex, environment, batch, and precomputed PCs.

    Returns DataFrame with marker index, N, beta, SE, t, p, R2 and MAF.
    Missing genotype rows are omitted marker-by-marker. Genotype mean
    imputation is deliberately not used in association testing.
    """
    G=np.asarray(genotypes,dtype=float); y=np.asarray(phenotype,dtype=float)
    if G.ndim!=2 or y.ndim!=1 or G.shape[0]!=len(y): raise ValueError("genotypes must be samples x markers and match phenotype")
    base=_design_matrix(covariates,len(y)); rows=[]
    for j in range(G.shape[1]):
        x=G[:,j]; keep=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(base),axis=1)
        xj=x[keep]; yj=y[keep]; C=base[keep]
        n=len(yj)
        if n<max(min_samples,C.shape[1]+3) or len(np.unique(xj))<2:
            rows.append({"marker_index":j,"n":n,"beta":np.nan,"se":np.nan,"t_stat":np.nan,"p_value":np.nan,"r_squared":np.nan,"maf":np.nan})
            continue
        X=np.column_stack([C,xj])
        rank=np.linalg.matrix_rank(X)
        if rank<X.shape[1]:
            rows.append({"marker_index":j,"n":n,"beta":np.nan,"se":np.nan,"t_stat":np.nan,"p_value":np.nan,"r_squared":np.nan,"maf":float(min(np.mean(xj)/2,1-np.mean(xj)/2))})
            continue
        coef,_,_,_=np.linalg.lstsq(X,yj,rcond=None)
        resid=yj-X@coef; df=n-rank
        s2=float(resid@resid/df)
        cov=s2*np.linalg.inv(X.T@X)
        se=float(np.sqrt(max(cov[-1,-1],0)))
        tstat=float(coef[-1]/se) if se>0 else np.nan
        pval=float(2*student_t.sf(abs(tstat),df)) if np.isfinite(tstat) else np.nan
        tss=float(np.sum((yj-yj.mean())**2))
        r2=float(1-(resid@resid)/tss) if tss>0 else np.nan
        p=float(np.mean(xj)/2)
        rows.append({"marker_index":j,"n":n,"beta":float(coef[-1]),"se":se,"t_stat":tstat,"p_value":pval,"r_squared":r2,"maf":min(p,1-p)})
    result=pd.DataFrame(rows)
    if len(result):
        result["p_adjusted_bonferroni"]=np.minimum(result["p_value"]*len(result),1.0)
        result["neg_log10_p"]=-np.log10(result["p_value"].clip(lower=np.finfo(float).tiny))
    return result

def genomic_inflation_factor(p_values):
    """Approximate lambda_GC from 1-df chi-square statistics derived from p-values."""
    from scipy.stats import chi2
    p=np.asarray(p_values,dtype=float); p=p[np.isfinite(p)&(p>0)&(p<=1)]
    if len(p)==0: return float("nan")
    chisq=chi2.isf(p,1)
    return float(np.median(chisq)/0.4549364)

def benjamini_hochberg(p_values):
    """Benjamini-Hochberg adjusted p-values, returned in original order."""
    p=np.asarray(p_values,dtype=float); out=np.full(p.shape,np.nan)
    idx=np.where(np.isfinite(p)&(p>=0)&(p<=1))[0]
    if len(idx)==0:return out
    order=idx[np.argsort(p[idx])]; m=len(order)
    q=p[order]*m/np.arange(1,m+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    out[order]=np.minimum(q,1)
    return out
