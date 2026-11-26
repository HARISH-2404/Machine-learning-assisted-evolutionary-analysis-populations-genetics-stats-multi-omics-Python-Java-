"""Gene diversity, effective allele number and multilocus summaries."""
import numpy as np
import pandas as pd
from .heterozygosity import heterozygosity
from ._core import save

def diversity(genotypes,population_col=None,id_col="genotype_id",output=None):
    h=heterozygosity(genotypes,population_col,id_col)
    h["effective_allele_number"]=1/(1-h.expected_heterozygosity).replace(0,np.nan)
    if output: save(h,output)
    return h

def diversity_summary(diversity_table):
    return diversity_table.groupby("population_id",as_index=False).agg(
        mean_Ho=("observed_heterozygosity","mean"),
        mean_He=("expected_heterozygosity","mean"),
        mean_effective_alleles=("effective_allele_number","mean"),
        loci=("marker_id","count"))
