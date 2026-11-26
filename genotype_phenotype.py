"""Genotype–phenotype integration with location/time/generation metadata."""
import pandas as pd

def integrate_genotype_phenotype(genotypes, phenotypes, metadata=None, join="inner"):
    """Join genotype and phenotype tables by sample ID.

    Phenotype rows may represent plot/replicate observations; provide a metadata
    table to retain location, season, time point and generation columns.
    """
    g=genotypes.copy(); p=phenotypes.copy()
    g.index=g.index.astype(str); p.index=p.index.astype(str)
    if metadata is not None:
        m=metadata.copy()
        if "sample_id" not in m.columns: raise ValueError("metadata needs sample_id")
        p=p.reset_index(names="sample_id").merge(m,on="sample_id",how="left").set_index("sample_id")
    if join not in {"inner","outer"}: raise ValueError("join must be inner or outer")
    out=g.join(p,how=join,lsuffix="__geno",rsuffix="__pheno")
    return out
