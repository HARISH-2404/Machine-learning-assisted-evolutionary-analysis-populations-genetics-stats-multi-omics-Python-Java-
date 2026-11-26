"""Random Forest wrappers."""
def fit_random_forest(X,y,task="regression",n_estimators=500,max_depth=None,seed=42,**kwargs):
    try:
        from sklearn.ensemble import RandomForestRegressor,RandomForestClassifier
    except ImportError as e: raise ImportError("pip install scikit-learn") from e
    if task=="regression": model=RandomForestRegressor(n_estimators=n_estimators,max_depth=max_depth,random_state=seed,n_jobs=-1,**kwargs)
    elif task=="classification": model=RandomForestClassifier(n_estimators=n_estimators,max_depth=max_depth,random_state=seed,n_jobs=-1,class_weight="balanced",**kwargs)
    else: raise ValueError("task must be regression or classification")
    model.fit(X,y); return {"model":model,"predictions":model.predict(X),"feature_importances":model.feature_importances_}
