"""Single-QTL interval scan and simple marker-trait regression.

The interval scan here is a marker-grid scan: each marker is tested as a
putative QTL position. It is not composite interval mapping and does not
calculate recombination-based LOD profiles between flanking markers.
"""
import numpy as np
import pandas as pd
from scipy.stats import f as f_dist

def marker_qtl_scan(genotypes, phenotype, marker_map, covariates=None, min_samples=10):
    """Perform marker-wise additive QTL scan; returns F statistic and LOD.

    genotypes shape is samples x markers (0/1/2); marker_map requires
    chromosome and position columns and has one row per marker.
    """
    G=np.asarray(genotypes,dtype=float); y=np.asarray(phenotype,dtype=float)
    mm=pd.DataFrame(marker_map).reset_index(drop=True)
    if G.ndim!=2 or G.shape[0]!=len(y) or G.shape[1]!=len(mm): raise ValueError("genotype dimensions must match phenotype and marker map")
    C=np.ones((len(y),1)) if covariates is None else np.column_stack([np.ones(len(y)),np.asarray(covariates,dtype=float)])
    if C.ndim!=2 or C.shape[0]!=len(y): raise ValueError("covariates rows mismatch")
    rows=[]
    for j in range(G.shape[1]):
        x=G[:,j]; keep=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(C),axis=1)
        X0=C[keep]; yy=y[keep]; xx=x[keep]; n=len(yy)
        if n<max(min_samples,X0.shape[1]+4) or len(np.unique(xx))<2:
            fval=lod=pval=np.nan
        else:
            X1=np.column_stack([X0,xx])
            if np.linalg.matrix_rank(X1)<X1.shape[1]: fval=lod=pval=np.nan
            else:
                b0=np.linalg.lstsq(X0,yy,rcond=None)[0]; rss0=float(np.sum((yy-X0@b0)**2))
                b1=np.linalg.lstsq(X1,yy,rcond=None)[0]; rss1=float(np.sum((yy-X1@b1)**2))
                df1=1; df2=n-np.linalg.matrix_rank(X1)
                fval=((rss0-rss1)/df1)/(rss1/df2) if rss1>0 else np.inf
                pval=float(f_dist.sf(fval,df1,df2))
                lod=(n/2)*np.log10(rss0/rss1) if rss1>0 and rss0>0 else np.nan
        row=mm.iloc[j].to_dict()
        row.update({"marker_index":j,"n":n,"F_statistic":fval,"LOD":lod,"p_value":pval})
        rows.append(row)
    return pd.DataFrame(rows)

def permutation_lod_threshold(genotypes, phenotype, marker_map, n_permutations=1000, alpha=0.05, seed=1, covariates=None):
    """Estimate genome-wide LOD threshold by phenotype permutation.

    For covariate-adjusted analyses this simple permutation is not a substitute
    for a carefully chosen residual/permutation scheme.
    """
    if n_permutations<1 or not 0<alpha<1: raise ValueError("invalid permutation settings")
    rng=np.random.default_rng(seed); maxima=[]
    y=np.asarray(phenotype,dtype=float)
    for _ in range(int(n_permutations)):
        scan=marker_qtl_scan(genotypes,rng.permutation(y),marker_map,covariates=covariates)
        vals=scan["LOD"].to_numpy(dtype=float); vals=vals[np.isfinite(vals)]
        maxima.append(float(vals.max()) if len(vals) else np.nan)
    a=np.asarray(maxima); a=a[np.isfinite(a)]
    return float(np.quantile(a,1-alpha)) if len(a) else float("nan")
