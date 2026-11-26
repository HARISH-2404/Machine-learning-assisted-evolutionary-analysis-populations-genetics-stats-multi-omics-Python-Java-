"""Transcriptome-wide association analysis (TWAS) utilities."""
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

def twas_expression_association(expression, phenotype, covariates=None, min_samples=10):
    """Associate each gene's expression with a quantitative phenotype.

    This is expression-trait regression, not a genetically predicted-expression
    TWAS. It does not establish causality and may be confounded by environment,
    tissue composition or reverse causation.
    """
    E=np.asarray(expression,float); y=np.asarray(phenotype,float)
    if E.ndim!=2 or y.ndim!=1 or E.shape[0]!=len(y):
        raise ValueError("expression must be samples x genes and match phenotype")
    C=np.ones((len(y),1)) if covariates is None else np.column_stack([np.ones(len(y)),np.asarray(covariates,float)])
    rows=[]
    for j in range(E.shape[1]):
        x=E[:,j]; keep=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(C),axis=1)
        xx=x[keep]; yy=y[keep]; cc=C[keep]; n=len(yy)
        if n<max(min_samples,cc.shape[1]+3) or len(np.unique(xx))<2:
            beta=se=tstat=pval=np.nan
        else:
            X=np.column_stack([cc,xx])
            if np.linalg.matrix_rank(X)<X.shape[1]: beta=se=tstat=pval=np.nan
            else:
                b=np.linalg.lstsq(X,yy,rcond=None)[0]; resid=yy-X@b; df=n-X.shape[1]
                cov=(resid@resid/df)*np.linalg.inv(X.T@X)
                beta=float(b[-1]); se=float(np.sqrt(max(cov[-1,-1],0)))
                tstat=beta/se if se>0 else np.nan
                pval=float(2*student_t.sf(abs(tstat),df)) if np.isfinite(tstat) else np.nan
        rows.append({"gene_id":f"gene_{j+1}","n":n,"beta":beta,"se":se,"t_statistic":tstat,"p_value":pval})
    out=pd.DataFrame(rows)
    if len(out):
        p=out.p_value.to_numpy(); ok=np.isfinite(p); out["p_adjusted_BH"]=np.nan
        if ok.any():
            ids=np.where(ok)[0]; order=ids[np.argsort(p[ok])]; m=len(order)
            q=np.minimum.accumulate((p[order]*m/np.arange(1,m+1))[::-1])[::-1]
            out.loc[order,"p_adjusted_BH"]=np.minimum(q,1)
    return out
