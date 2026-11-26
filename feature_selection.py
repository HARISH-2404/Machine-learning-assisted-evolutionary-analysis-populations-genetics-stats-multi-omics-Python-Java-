"""Simple univariate feature scoring and top-k selection."""
import numpy as np

def univariate_feature_scores(X,y):
    X=np.asarray(X,float); y=np.asarray(y,float)
    if X.ndim!=2 or y.ndim!=1 or len(y)!=len(X): raise ValueError("Invalid dimensions")
    yc=y-y.mean(); den_y=np.sqrt((yc*yc).sum()); scores=[]
    for j in range(X.shape[1]):
        x=X[:,j]-X[:,j].mean(); den=np.sqrt((x*x).sum()* (den_y**2))
        scores.append(abs(float(x@yc/den)) if den else 0.0)
    return np.asarray(scores)

def select_top_features(X,y,k):
    scores=univariate_feature_scores(X,y)
    if k<1 or k>len(scores): raise ValueError("k must be between 1 and number of features")
    idx=np.argsort(scores)[::-1][:k]
    return {"indices":idx,"scores":scores[idx],"X_selected":np.asarray(X)[:,idx]}
