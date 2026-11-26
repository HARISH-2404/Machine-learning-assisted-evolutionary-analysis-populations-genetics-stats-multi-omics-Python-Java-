"""Multi-trait ridge prediction for several correlated phenotypes.

Fits one ridge model per trait using a shared marker matrix. This is a
multi-output baseline; it does not estimate residual/genetic covariance
matrices or implement a full multivariate mixed model.
"""
import numpy as np

def multi_trait_rrblup(X, Y, lambda_=1.0, standardize=True):
    X=np.asarray(X,dtype=float); Y=np.asarray(Y,dtype=float)
    if X.ndim!=2 or Y.ndim!=2 or X.shape[0]!=Y.shape[0] or lambda_<=0:
        raise ValueError("X must be n x p, Y n x traits, and lambda positive")
    if not np.isfinite(X).all() or not np.isfinite(Y).all(): raise ValueError("Inputs must be finite")
    xm=X.mean(axis=0); xs=X.std(axis=0); xs[xs==0]=1
    Z=(X-xm)/xs if standardize else X.copy()
    ym=Y.mean(axis=0); Yc=Y-ym
    B=np.linalg.solve(Z.T@Z+lambda_*np.eye(Z.shape[1]),Z.T@Yc)
    pred=ym+Z@B
    return {"marker_effects":B,"trait_intercepts":ym,"predictions":pred,
            "observed_trait_correlation":np.corrcoef(Y,rowvar=False)}
