"""Regression evaluation metrics."""
import numpy as np

def regression_metrics(y_true,y_pred):
    y=np.asarray(y_true,float); p=np.asarray(y_pred,float)
    if y.shape!=p.shape or y.ndim!=1 or len(y)==0: raise ValueError("Equal nonempty vectors required")
    err=y-p; mse=float(np.mean(err**2)); den=float(np.sum((y-y.mean())**2))
    return {"mse":mse,"rmse":mse**0.5,"mae":float(np.mean(np.abs(err))),
            "r_squared":1-float(np.sum(err**2))/den if den else float("nan")}
