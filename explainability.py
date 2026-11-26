"""Model-agnostic permutation importance on a held-out dataset."""
import numpy as np
from .model_evaluation import regression_metrics

def permutation_importance(predict_fn,X,y,n_repeats=10,seed=42):
    X=np.asarray(X,float); y=np.asarray(y,float)
    baseline=regression_metrics(y,predict_fn(X))["rmse"]
    rng=np.random.default_rng(seed); importances=[]
    for j in range(X.shape[1]):
        scores=[]
        for _ in range(n_repeats):
            xp=X.copy(); xp[:,j]=rng.permutation(xp[:,j])
            scores.append(regression_metrics(y,predict_fn(xp))["rmse"]-baseline)
        importances.append({"feature":j,"importance_mean":float(np.mean(scores)),
                            "importance_sd":float(np.std(scores,ddof=1)) if n_repeats>1 else 0.0})
    return sorted(importances,key=lambda r:r["importance_mean"],reverse=True)
