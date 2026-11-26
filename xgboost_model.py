"""XGBoost wrapper."""
def fit_xgboost(X,y,task="regression",n_estimators=500,max_depth=5,learning_rate=.03,seed=42,**kwargs):
    try: from xgboost import XGBRegressor,XGBClassifier
    except ImportError as e: raise ImportError("pip install xgboost") from e
    common=dict(n_estimators=n_estimators,max_depth=max_depth,learning_rate=learning_rate,random_state=seed,n_jobs=-1,**kwargs)
    if task=="regression": model=XGBRegressor(objective="reg:squarederror",**common)
    elif task=="classification": model=XGBClassifier(objective="binary:logistic",eval_metric="logloss",**common)
    else: raise ValueError("task must be regression or classification")
    model.fit(X,y); return {"model":model,"predictions":model.predict(X),"feature_importances":model.feature_importances_}
