"""VCF parsing, genotype dosage matrix and variant-level summaries."""
import math
import numpy as np
import pandas as pd
from ._core import parse_vcf, parse_genotypes

def load_vcf_dosage(path, snps_only=False, biallelic_only=True):
    """Return (variant_metadata, dosage_dataframe); rows are variants, columns samples.

    Dosage is alternate allele count 0/1/2. This parser is for diploid biallelic
    records; it does not infer polyploid dosage or phase.
    """
    metadata=[]; rows=[]; samples=None
    for r in parse_vcf(path):
        if snps_only and (len(r["ref"])!=1 or len(r["alt"])!=1): continue
        if biallelic_only and "," in r["alt"]: continue
        d=parse_genotypes(r)
        if samples is None: samples=r["sample_names"]
        metadata.append({k:r[k] for k in ("chrom","pos","id","ref","alt","qual","filter")})
        rows.append(d)
    if samples is None: samples=[]
    ids=[m["id"] if m["id"] not in ("",".") else f'{m["chrom"]}:{m["pos"]}' for m in metadata]
    return pd.DataFrame(metadata), pd.DataFrame(rows,index=ids,columns=samples,dtype=float)

def variant_statistics(dosage):
    """Per-variant call rate, alternate frequency, MAF and observed heterozygosity."""
    x=np.asarray(dosage,dtype=float)
    if x.ndim!=2: raise ValueError("dosage must be variants x samples")
    out=[]
    for i,row in enumerate(x):
        called=np.isfinite(row); n=int(called.sum()); nmiss=len(row)-n
        p=float(np.nansum(row)/(2*n)) if n else float("nan")
        maf=min(p,1-p) if n else float("nan")
        ho=float(np.mean(row[called]==1)) if n else float("nan")
        out.append({"variant_index":i,"n_called":n,"n_missing":nmiss,
                    "call_rate":n/len(row) if len(row) else float("nan"),
                    "alt_frequency":p,"maf":maf,"observed_heterozygosity":ho})
    return out

def filter_variants(dosage, min_call_rate=0.9, min_maf=0.01, max_maf=0.5):
    """Return boolean keep mask; thresholds are inclusive."""
    if not 0<=min_call_rate<=1 or not 0<=min_maf<=max_maf<=0.5:
        raise ValueError("Invalid filtering thresholds")
    stats=variant_statistics(dosage)
    return np.array([s["call_rate"]>=min_call_rate and min_maf<=s["maf"]<=max_maf for s in stats])
