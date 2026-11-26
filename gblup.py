"""Genomic BLUP using a supplied genomic relationship matrix (GRM).

Solves y = 1*mu + u + e, u ~ N(0, G*sigma_g2), e ~ N(0, I*sigma_e2).
For fixed variance ratio lambda=sigma_e2/sigma_g2, the returned solution is
a penalized kernel predictor. Variance components are not estimated by REML.
"""
import numpy as np

def gblup(G, y, lambda_=1.0, fixed_effects=None):
    G=np.asarray(G,dtype=float); y=np.asarray(y,dtype=float)
    n=len(y)
    if G.shape!=(n,n) or not np.isfinite(G).all() or not np.isfinite(y).all():
        raise ValueError("G must be finite n x n and y a finite length-n vector")
    if lambda_<=0: raise ValueError("lambda_ must be positive")
    if not np.allclose(G,G.T,atol=1e-8): raise ValueError("G must be symmetric")
    X=np.ones((n,1)) if fixed_effects is None else np.column_stack([np.ones(n),np.asarray(fixed_effects,dtype=float)])
    if X.shape[0]!=n: raise ValueError("fixed_effects row count must match y")
    # Mixed-model equations in kernel form, with fixed effects profiled by GLS.
    V=G+lambda_*np.eye(n)
    Vi=np.linalg.solve(V,np.eye(n))
    XtViX=X.T@Vi@X
    beta=np.linalg.solve(XtViX,X.T@Vi@y)
    u=G@Vi@(y-X@beta)
    pred=X@beta+u
    return {"fixed_effects":beta,"genomic_values":u,"predictions":pred,
            "residuals":y-pred,"lambda":float(lambda_)}
