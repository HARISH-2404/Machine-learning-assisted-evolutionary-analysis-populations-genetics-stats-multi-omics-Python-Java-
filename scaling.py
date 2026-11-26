"""Feature-block scaling; fit on training data only to prevent leakage."""
import numpy as np
import pandas as pd

def scale_feature_block(train, test=None, method="zscore"):
    """Scale columns using training statistics; returns (train_scaled, test_scaled)."""
    if method not in {"zscore", "minmax", "none"}: raise ValueError("Unsupported scaling method")
    a=train.astype(float).copy()
    b=None if test is None else test.astype(float).copy()
    if method == "none": return a,b
    if method == "zscore":
        center=a.mean(axis=0); spread=a.std(axis=0, ddof=0).replace(0,1)
        a=(a-center)/spread
        if b is not None: b=(b-center)/spread
    else:
        lo=a.min(axis=0); span=(a.max(axis=0)-lo).replace(0,1)
        a=(a-lo)/span
        if b is not None: b=(b-lo)/span
    return a,b
