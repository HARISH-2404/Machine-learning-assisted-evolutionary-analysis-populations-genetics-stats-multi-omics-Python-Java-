"""Join genotype and expression features at matched sample/gene level."""
import pandas as pd

def integrate_genotype_expression(genotypes, expression, sample_map=None, gene_map=None, join="inner"):
    """Align sample rows and optionally map variant IDs to genes.

    Inputs are sample-by-feature DataFrames. `gene_map` maps genotype feature IDs
    to gene IDs; unmapped variants retain their original IDs.
    """
    g=genotypes.copy(); e=expression.copy()
    if sample_map is not None:
        g=g.rename(index=sample_map); e=e.rename(index=sample_map)
    if gene_map is not None:
        g=g.rename(columns={k:v for k,v in gene_map.items()})
        # Multiple variants may map to one gene; aggregate dosage by mean.
        g=g.T.groupby(level=0).mean().T
    shared=g.index.intersection(e.index)
    if join=="inner":
        g=g.loc[shared]; e=e.loc[shared]
    elif join=="outer":
        g,e=g.align(e,join="outer",axis=0)
    else: raise ValueError("join must be inner or outer")
    g.columns=[f"geno::{c}" for c in g.columns]
    e.columns=[f"expr::{c}" for c in e.columns]
    return pd.concat([g,e],axis=1)
