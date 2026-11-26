"""Gene co-expression correlation networks."""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

def coexpression_matrix(expression, method="pearson", min_periods=4):
    """Return gene x gene correlation matrix from samples x genes expression."""
    df=expression.copy() if isinstance(expression,pd.DataFrame) else pd.DataFrame(expression)
    if method not in ("pearson","spearman"): raise ValueError("method must be pearson or spearman")
    return df.corr(method=method,min_periods=min_periods)

def coexpression_edges(expression, threshold=0.7, method="pearson", absolute=True):
    """Return upper-triangle gene pairs exceeding correlation threshold."""
    corr=coexpression_matrix(expression,method=method)
    genes=list(corr.columns); rows=[]
    for i in range(len(genes)):
        for j in range(i+1,len(genes)):
            r=corr.iloc[i,j]
            if np.isfinite(r) and (abs(r)>=threshold if absolute else r>=threshold):
                rows.append({"gene_1":genes[i],"gene_2":genes[j],"correlation":float(r)})
    return pd.DataFrame(rows,columns=["gene_1","gene_2","correlation"])

def module_membership(expression, module_gene_lists):
    """Compute mean within-module expression profile for each supplied gene list."""
    df=expression.copy() if isinstance(expression,pd.DataFrame) else pd.DataFrame(expression)
    out={}
    for name,genes in module_gene_lists.items():
        present=[g for g in genes if g in df.columns]
        if not present: continue
        out[name]=df[present].mean(axis=1)
    return pd.DataFrame(out,index=df.index)
