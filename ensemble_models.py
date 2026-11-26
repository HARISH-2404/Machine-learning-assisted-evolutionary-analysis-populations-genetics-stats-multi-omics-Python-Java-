"""ExtraTrees, gradient boosting, histogram gradient boosting, KNN and Naive Bayes."""
def fit_model(name,X,y,task="regression",seed=42,**kwargs):
    try:
        from sklearn.ensemble import ExtraTreesRegressor,ExtraTreesClassifier,GradientBoostingRegressor,GradientBoostingClassifier,HistGradientBoostingRegressor,HistGradientBoostingClassifier
        from sklearn.neighbors import KNeighborsRegressor,KNeighborsClassifier
        from sklearn.naive_bayes import GaussianNB
    except ImportError as e: raise ImportError("pip install scikit-learn") from e
    choices={"extra_trees":(ExtraTreesRegressor,ExtraTreesClassifier),"gradient_boosting":(GradientBoostingRegressor,GradientBoostingClassifier),"hist_gradient_boosting":(HistGradientBoostingRegressor,HistGradientBoostingClassifier),"knn":(KNeighborsRegressor,KNeighborsClassifier),"naive_bayes":(None,GaussianNB)}
    if name not in choices: raise ValueError("Unknown model")
    cls=choices[name][0 if task=="regression" else 1]
    if cls is None: raise ValueError("Naive Bayes is classification-only")
    model=cls(random_state=seed,**kwargs) if name in ("extra_trees","gradient_boosting") else cls(**kwargs)
    model.fit(X,y); out={"model":model,"predictions":model.predict(X)}
    if hasattr(model,"feature_importances_"): out["feature_importances"]=model.feature_importances_
    return out
