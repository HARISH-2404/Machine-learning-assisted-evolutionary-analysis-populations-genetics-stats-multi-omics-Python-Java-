"""RR-BLUP marker-effect prediction."""
import numpy as np

def rrblup(X, y, lambda_=1.0, standardize=True):
    """Estimate marker effects from (X'X + lambda I)b = X'y."""
    X=np.asarray(X,dtype=float); y=np.asarray(y,dtype=float)
    if X.ndim!=2 or y.ndim!=1 or X.shape[0]!=len(y) or lambda_<=0:
        raise ValueError("Invalid X/y dimensions or lambda")
    if not np.isfinite(X).all() or not np.isfinite(y).all(): raise ValueError("Inputs must be finite")
    means=X.mean(axis=0); scales=X.std(axis=0); scales[scales==0]=1
    Z=(X-means)/scales if standardize else X.copy()
    ym=y.mean(); yc=y-ym
    effects=np.linalg.solve(Z.T@Z+lambda_*np.eye(Z.shape[1]),Z.T@yc)
    fitted=ym+Z@effects
    return {"marker_effects":effects,"intercept":ym,"predictions":fitted,
            "marker_means":means,"marker_scales":scales}
