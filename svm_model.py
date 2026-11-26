"""SVM regression/classification wrapper with standardization."""
def fit_svm(X,y,task="regression",kernel="rbf",C=1.0,epsilon=.1,gamma="scale"):
    try:
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVR,SVC
    except ImportError as e: raise ImportError("pip install scikit-learn") from e
    est=SVR(kernel=kernel,C=C,epsilon=epsilon,gamma=gamma) if task=="regression" else SVC(kernel=kernel,C=C,gamma=gamma,probability=True,class_weight="balanced")
    model=make_pipeline(StandardScaler(),est); model.fit(X,y); return {"model":model,"predictions":model.predict(X)}
