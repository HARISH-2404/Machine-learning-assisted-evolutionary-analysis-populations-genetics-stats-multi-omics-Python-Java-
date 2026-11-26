"""Expression-QTL association by per-marker linear regression.

This is a transparent single-marker model. For population structure, relatedness,
repeated measures and multiple environments, use a validated mixed-model eQTL
workflow (e.g. tensorQTL/FastQTL with appropriate covariates).
"""
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

def eqtl_scan(genotypes, expression, covariates=None, min_samples=10):
    """Test each genotype marker against each expression feature.

    genotypes: samples x markers (diploid dosage 0/1/2)
    expression: samples x genes (normalized/log expression)
    Returns long table with marker_id, gene_id, beta, SE, t, p and N.
    """
    G=np.asarray(genotypes,float); E=np.asarray(expression,float)
    if G.ndim!=2 or E.ndim!=2 or G.shape[0]!=E.shape[0]:
        raise ValueError("genotypes and expression must be 2D with matching sample rows")
    n=G.shape[0]
    C=np.ones((n,1)) if covariates is None else np.column_stack([np.ones(n),np.asarray(covariates,float)])
    if C.shape[0]!=n: raise ValueError("covariate rows mismatch")
    rows=[]
    for j in range(G.shape[1]):
        for k in range(E.shape[1]):
            x=G[:,j]; y=E[:,k]
            keep=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(C),axis=1)
            xx=x[keep]; yy=y[keep]; cc=C[keep]; N=len(yy)
            if N<max(min_samples,cc.shape[1]+3) or len(np.unique(xx))<2:
                beta=se=tstat=pval=np.nan
            else:
                X=np.column_stack([cc,xx])
                if np.linalg.matrix_rank(X)<X.shape[1]: beta=se=tstat=pval=np.nan
                else:
                    b=np.linalg.lstsq(X,yy,rcond=None)[0]; r=yy-X@b; df=N-X.shape[1]
                    s2=float(r@r/df); cov=s2*np.linalg.inv(X.T@X)
                    beta=float(b[-1]); se=float(np.sqrt(max(cov[-1,-1],0)))
                    tstat=beta/se if se>0 else np.nan
                    pval=float(2*student_t.sf(abs(tstat),df)) if np.isfinite(tstat) else np.nan
            rows.append({"marker_id":f"marker_{j+1}","gene_id":f"gene_{k+1}","n":N,
                         "beta":beta,"se":se,"t_statistic":tstat,"p_value":pval})
    out=pd.DataFrame(rows)
    if len(out):
        out["p_adjusted_BH"]=np.nan
        for gene,idx in out.groupby("gene_id").groups.items():
            p=out.loc[idx,"p_value"].to_numpy(); valid=np.isfinite(p)
            if valid.any():
                order=np.argsort(p[valid]); pv=p[valid][order]; m=len(pv)
                q=np.minimum.accumulate((pv*m/np.arange(1,m+1))[::-1])[::-1]
                loc=np.asarray(list(idx))[valid][order]; out.loc[loc,"p_adjusted_BH"]=np.minimum(q,1)
    return out
