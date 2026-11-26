"""Ridge-regression BLUP approximation for marker effects.

For centered/scaled marker matrix X and phenotype y, solves
b=(X'X + lambda I)^-1 X'y. This is equivalent to ridge regression and can
approximate genomic prediction under specific assumptions; it is not a full
multi-environment mixed model or REML variance-component estimator.
"""
import numpy as np

def ridge_blup(X, y, alpha=1.0, standardize=True):
    X=np.asarray(X,dtype=float); y=np.asarray(y,dtype=float)
    if X.ndim!=2 or y.ndim!=1 or X.shape[0]!=len(y) or alpha<=0:
        raise ValueError("Invalid X/y dimensions or alpha")
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("X and y must be finite")
    means=X.mean(axis=0); scales=X.std(axis=0,ddof=0)
    scales[scales==0]=1.0
    Z=(X-means)/scales if standardize else X.copy()
    ym=y.mean(); yc=y-ym
    beta=np.linalg.solve(Z.T@Z+alpha*np.eye(Z.shape[1]),Z.T@yc)
    fitted=ym+Z@beta
    return {"marker_effects":beta,"intercept":ym,"predictions":fitted,
            "r_squared":1-float(np.sum((y-fitted)**2))/float(np.sum((y-ym)**2)) if np.sum((y-ym)**2)>0 else float("nan"),
            "marker_means":means,"marker_scales":scales}
