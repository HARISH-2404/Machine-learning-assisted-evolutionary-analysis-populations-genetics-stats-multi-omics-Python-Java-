"""Within-population inbreeding coefficient FIS = 1 - Ho/He.

Reports per-population multilocus ratio of sums, avoiding an unweighted
average of unstable locus-wise ratios. Also returns locus-level values.
"""
import numpy as np
import pandas as pd
from .heterozygosity import heterozygosity
from ._core import save

def fis(genotypes,population_col="population_id",id_col="genotype_id",output=None):
    h=heterozygosity(genotypes,population_col,id_col)
    h["FIS_locus"]=1-h.observed_heterozygosity/h.expected_heterozygosity.replace(0,np.nan)
    rows=[]
    for pop,g in h.groupby("population_id"):
        ho=g.observed_heterozygosity.sum(min_count=1)
        he=g.expected_heterozygosity.sum(min_count=1)
        rows.append({"population_id":pop,"Ho_sum":ho,"He_sum":he,
                     "FIS_ratio_of_sums":1-ho/he if pd.notna(he) and he>0 else np.nan,
                     "loci":len(g)})
    summary=pd.DataFrame(rows)
    if output:
        save(h, str(output).replace(".csv","_locus.csv"))
        save(summary,output)
    return h,summary
