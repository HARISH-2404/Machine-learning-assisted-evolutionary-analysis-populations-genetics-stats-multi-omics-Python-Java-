"""Linear regression baselines with ridge and coordinate-descent lasso."""
import numpy as np

def ridge_regression(X,y,alpha=1.0,standardize=True):
    X=np.asarray(X,float); y=np.asarray(y,float)
    if X.ndim!=2 or y.ndim!=1 or len(y)!=X.shape[0] or alpha<0: raise ValueError("Invalid inputs")
    if not np.isfinite(X).all() or not np.isfinite(y).all(): raise ValueError("Inputs must be finite")
    mu=X.mean(0); sd=X.std(0); sd[sd==0]=1
    Z=(X-mu)/sd if standardize else X.copy(); ym=y.mean(); yc=y-ym
    b=np.linalg.solve(Z.T@Z+alpha*np.eye(Z.shape[1]),Z.T@yc)
    pred=ym+Z@b
    return {"coefficients":b,"intercept":ym,"predictions":pred,"feature_mean":mu,"feature_scale":sd}

def lasso_coordinate_descent(X,y,alpha=0.1,max_iter=5000,tol=1e-7):
    X=np.asarray(X,float); y=np.asarray(y,float)
    if X.ndim!=2 or y.ndim!=1 or len(y)!=X.shape[0] or alpha<0: raise ValueError("Invalid inputs")
    mu=X.mean(0); sd=X.std(0); sd[sd==0]=1; Z=(X-mu)/sd; ym=y.mean(); yc=y-ym
    n,p=Z.shape; b=np.zeros(p); norms=(Z*Z).sum(0)/n
    for _ in range(max_iter):
        old=b.copy()
        for j in range(p):
            residual=yc-Z@b+Z[:,j]*b[j]
            rho=(Z[:,j]@residual)/n
            if norms[j]==0: b[j]=0
            elif rho>alpha: b[j]=(rho-alpha)/norms[j]
            elif rho < -alpha: b[j]=(rho+alpha)/norms[j]
            else: b[j]=0
        if np.max(np.abs(b-old))<tol: break
    pred=ym+Z@b
    return {"coefficients":b,"intercept":ym,"predictions":pred,"iterations":_+1}
