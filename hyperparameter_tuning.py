"""Small transparent grid search for ridge regression."""
import numpy as np
from .cross_validation import kfold_indices
from .regression import ridge_regression

def grid_search_ridge(X,y,alphas=(0.01,0.1,1,10,100),n_splits=5,groups=None,seed=42):
    X=np.asarray(X,float); y=np.asarray(y,float)
    if groups is None: folds=list(kfold_indices(len(y),n_splits,True,seed))
    else:
        from .cross_validation import grouped_kfold_indices
        folds=list(grouped_kfold_indices(groups,n_splits,seed))
    results=[]
    for a in alphas:
        scores=[]
        for tr,te in folds:
            # Fit preprocessing only on training data to avoid leakage.
            fit=ridge_regression(X[tr],y[tr],alpha=float(a))
            mu=fit["feature_mean"]; sd=fit["feature_scale"]
            Z=(X[te]-mu)/sd
            pred=fit["intercept"]+Z@fit["coefficients"]
            scores.append(float(np.mean((y[te]-pred)**2)))
        results.append({"alpha":float(a),"mean_cv_mse":float(np.mean(scores))})
    best=min(results,key=lambda r:r["mean_cv_mse"])
    return {"best_alpha":best["alpha"],"best_cv_mse":best["mean_cv_mse"],"results":results}
