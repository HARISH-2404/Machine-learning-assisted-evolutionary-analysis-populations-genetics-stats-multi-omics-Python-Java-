"""Population differentiation and diversity estimators for biallelic loci."""
import numpy as np
import pandas as pd

def expected_heterozygosity(p):
    p=np.asarray(p,dtype=float)
    return 2*p*(1-p)

def pairwise_fst_from_frequencies(freq_table, population_col="population", locus_cols=("chromosome","position","gene_id"),
                                  frequency_col="allele_frequency"):
    """Simple per-locus heterozygosity-based FST approximation for two populations.

    FST=(HT-HS)/HT, with HT based on mean allele frequency and HS on mean
    within-population heterozygosity. Use sample-size-corrected estimators
    (e.g. Weir-Cockerham) for publication-grade inference.
    """
    pops=list(freq_table[population_col].dropna().unique())
    rows=[]
    for i in range(len(pops)):
        for j in range(i+1,len(pops)):
            a=freq_table[freq_table[population_col]==pops[i]]
            b=freq_table[freq_table[population_col]==pops[j]]
            m=a.merge(b,on=list(locus_cols),suffixes=("_a","_b"))
            for _,r in m.iterrows():
                p1=float(r[frequency_col+"_a"]); p2=float(r[frequency_col+"_b"])
                ht=2*((p1+p2)/2)*(1-(p1+p2)/2)
                hs=(2*p1*(1-p1)+2*p2*(1-p2))/2
                fst=(ht-hs)/ht if ht>0 else np.nan
                rows.append({"population_1":pops[i],"population_2":pops[j],
                             **{c:r[c] for c in locus_cols},"fst_approx":fst})
    return pd.DataFrame(rows)

def summarize_diversity(freq_table, group_cols=("population","location"), frequency_col="allele_frequency"):
    x=freq_table.copy()
    x["heterozygosity"]=expected_heterozygosity(x[frequency_col])
    return x.groupby(list(group_cols),dropna=False)["heterozygosity"].mean().reset_index(name="mean_expected_heterozygosity")
