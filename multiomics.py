"""Named multi-omics blocks and integration summaries."""
import pandas as pd
from .merge import merge_omics_blocks
from .association import pairwise_feature_correlations

def integrate_multiomics(blocks, join="inner"):
    """Return merged matrix plus per-block sample/feature coverage report."""
    merged=merge_omics_blocks(blocks,join=join,prefixes=True)
    report=pd.DataFrame([{"omics":name,"n_samples":len(df.index.unique()),
                          "n_features":df.shape[1],
                          "missing_fraction":float(df.isna().to_numpy().mean())}
                         for name,df in blocks.items()])
    return merged,report

def cross_block_associations(blocks, method="spearman", min_pairs=5):
    """Calculate pairwise feature correlations between every pair of omics blocks."""
    names=list(blocks); outputs=[]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            x=pairwise_feature_correlations(blocks[a],blocks[b],method,min_pairs)
            if not x.empty:
                x.insert(0,"block_b",b); x.insert(0,"block_a",a); outputs.append(x)
    return pd.concat(outputs,ignore_index=True) if outputs else pd.DataFrame(
        columns=["block_a","block_b","feature_a","feature_b","n","correlation"])
