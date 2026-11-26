"""Bayesian ridge/Gaussian marker regression using Gibbs sampling.

This module implements Bayesian ridge regression with Gaussian marker effects
and inverse-gamma priors on residual and marker-effect variances. It is a
transparent Bayesian baseline, not full BayesA/B/C/pi.
"""
import numpy as np

def bayesian_ridge_gibbs(X, y, iterations=3000, burn_in=1000, thin=5,
                         beta_var_prior=1.0, residual_var_prior=1.0,
                         prior_df=3.0, seed=123):
    X=np.asarray(X,dtype=float); y=np.asarray(y,dtype=float)
    if X.ndim!=2 or y.ndim!=1 or X.shape[0]!=len(y): raise ValueError("Invalid X/y dimensions")
    if iterations<=burn_in or burn_in<0 or thin<1: raise ValueError("Require iterations > burn_in >= 0 and thin >= 1")
    if min(beta_var_prior,residual_var_prior,prior_df)<=0: raise ValueError("Prior parameters must be positive")
    if not np.isfinite(X).all() or not np.isfinite(y).all(): raise ValueError("Inputs must be finite")
    n,p=X.shape; rng=np.random.default_rng(seed)
    xm=X.mean(axis=0); xs=X.std(axis=0); xs[xs==0]=1
    Z=(X-xm)/xs; ym=y.mean(); yc=y-ym
    b=np.zeros(p); sigma_e=max(float(np.var(yc)),1e-6); sigma_b=1.0
    xtx=Z.T@Z; xty=Z.T@yc
    samples_b=[]; samples_e=[]; samples_bvar=[]
    for it in range(iterations):
        precision=xtx/sigma_e+np.eye(p)/sigma_b
        cov=np.linalg.inv(precision)
        mean=cov@(xty/sigma_e)
        b=rng.multivariate_normal(mean,cov)
        resid=yc-Z@b
        # inverse-gamma draw via reciprocal gamma
        shape_e=(prior_df+n)/2
        scale_e=(residual_var_prior+resid@resid)/2
        sigma_e=1/rng.gamma(shape_e,1/scale_e)
        shape_b=(prior_df+p)/2
        scale_b=(beta_var_prior+b@b)/2
        sigma_b=1/rng.gamma(shape_b,1/scale_b)
        if it>=burn_in and (it-burn_in)%thin==0:
            samples_b.append(b.copy()); samples_e.append(sigma_e); samples_bvar.append(sigma_b)
    B=np.asarray(samples_b)
    if B.size==0: raise ValueError("No posterior samples retained")
    pred=ym+Z@B.mean(axis=0)
    return {"marker_effect_mean":B.mean(axis=0),"marker_effect_sd":B.std(axis=0,ddof=1) if len(B)>1 else np.zeros(p),
            "residual_variance_mean":float(np.mean(samples_e)),
            "marker_variance_mean":float(np.mean(samples_bvar)),
            "predictions":pred,"posterior_samples":len(B),
            "note":"Posterior marker variance is sampled internally; use posterior draws for full variance summaries."}
