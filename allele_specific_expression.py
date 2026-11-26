"""Allele-specific expression (ASE) count tests at heterozygous sites."""
import numpy as np
import pandas as pd
from scipy.stats import binomtest

def ase_binomial_test(counts, min_total=10, null_fraction=0.5):
    """Test reference-vs-alternate read balance using exact binomial tests.

    counts requires columns ref_count and alt_count. This basic test assumes
    independent reads and a 1:1 null. Mapping bias, overdispersion, copy number,
    RNA editing and allele-specific mapping must be addressed in real analyses.
    """
    df=counts.copy() if isinstance(counts,pd.DataFrame) else pd.DataFrame(counts)
    required={"ref_count","alt_count"}
    if not required.issubset(df.columns): raise ValueError("counts must contain ref_count and alt_count")
    rows=[]
    for idx,r in df.iterrows():
        ref=int(r.ref_count); alt=int(r.alt_count); total=ref+alt
        if ref<0 or alt<0: raise ValueError("read counts must be non-negative")
        if total<min_total:
            p=np.nan
        else:
            p=float(binomtest(alt,total,p=null_fraction,alternative="two-sided").pvalue)
        rows.append({"site_id":r.get("site_id",idx),"ref_count":ref,"alt_count":alt,
                     "total_count":total,"alt_fraction":alt/total if total else np.nan,"p_value":p})
    out=pd.DataFrame(rows)
    p=out.p_value.to_numpy(float); ok=np.isfinite(p); out["p_adjusted_BH"]=np.nan
    if ok.any():
        ids=np.where(ok)[0]; order=ids[np.argsort(p[ok])]; m=len(order)
        q=np.minimum.accumulate((p[order]*m/np.arange(1,m+1))[::-1])[::-1]
        out.loc[order,"p_adjusted_BH"]=np.minimum(q,1)
    return out
