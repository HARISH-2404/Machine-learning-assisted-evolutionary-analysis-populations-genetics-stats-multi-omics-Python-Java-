"""Cross-validation and nested hyperparameter validation."""
import numpy as np
def cross_validate_estimator(estimator,X,y,task="regression",cv=5,groups=None,scoring=None):
    try:
        from sklearn.model_selection import cross_validate,KFold,StratifiedKFold,GroupKFold
    except ImportError as e: raise ImportError("pip install scikit-learn") from e
    splitter=GroupKFold(cv) if groups is not None else (KFold(cv,shuffle=True,random_state=42) if task=="regression" else StratifiedKFold(cv,shuffle=True,random_state=42))
    result=cross_validate(estimator,X,y,groups=groups,cv=splitter,scoring=scoring,return_train_score=True)
    return {k:np.asarray(v).tolist() for k,v in result.items() if k.startswith(("test_","train_"))}
