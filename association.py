"""Cross-omics feature association screening."""
import numpy as np
import pandas as pd

def pairwise_feature_correlations(block_a, block_b, method="spearman", min_pairs=5):
    """Return long table of pairwise correlations for shared samples (screening only)."""
    common=block_a.index.intersection(block_b.index)
    a=block_a.loc[common].astype(float); b=block_b.loc[common].astype(float)
    rows=[]
    for ca in a.columns:
        for cb in b.columns:
            pair=pd.concat([a[ca],b[cb]],axis=1).dropna()
            if len(pair)<min_pairs: continue
            r=pair.iloc[:,0].corr(pair.iloc[:,1],method=method)
            rows.append((ca,cb,len(pair),r))
    return pd.DataFrame(rows,columns=["feature_a","feature_b","n","correlation"])
