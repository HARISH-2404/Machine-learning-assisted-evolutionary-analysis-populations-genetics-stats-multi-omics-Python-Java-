"""Population-genomic summaries commonly used in selective-sweep screening."""
import numpy as np
import pandas as pd

def site_frequency_spectrum(genotypes):
    """Folded SFS counts for diploid 0/1/2 dosage matrix (samples x variants)."""
    G=np.asarray(genotypes,dtype=float)
    if G.ndim!=2: raise ValueError("genotypes must be samples x variants")
    bins={}
    for j in range(G.shape[1]):
        x=G[:,j]; x=x[np.isfinite(x)]
        if len(x)==0: continue
        ac=int(np.sum(x)); an=2*len(x); minor=min(ac,an-ac)
        bins[minor]=bins.get(minor,0)+1
    return pd.DataFrame([{"minor_allele_count":k,"site_count":v} for k,v in sorted(bins.items())])

def windowed_heterozygosity(genotypes, positions, chromosomes=None, window_size=100000):
    """Calculate mean expected heterozygosity 2p(1-p) per genomic window."""
    G=np.asarray(genotypes,dtype=float); pos=np.asarray(positions,dtype=int)
    if G.ndim!=2 or G.shape[1]!=len(pos): raise ValueError("positions must match variants")
    chrvals=np.array(["1"]*len(pos) if chromosomes is None else chromosomes,dtype=str)
    if len(chrvals)!=len(pos): raise ValueError("chromosomes length mismatch")
    out=[]
    for chrom in np.unique(chrvals):
        js=np.where(chrvals==chrom)[0]
        if len(js)==0:continue
        lo=(int(pos[js].min())-1)//window_size*window_size+1
        hi=int(pos[js].max())
        for start in range(lo,hi+1,window_size):
            end=start+window_size-1; vals=[]
            for j in js[(pos[js]>=start)&(pos[js]<=end)]:
                x=G[:,j]; x=x[np.isfinite(x)]
                if len(x):
                    p=x.sum()/(2*len(x)); vals.append(2*p*(1-p))
            out.append({"chrom":chrom,"start":start,"end":end,"n_variants":len(vals),
                        "mean_heterozygosity":float(np.mean(vals)) if vals else np.nan})
    return pd.DataFrame(out)

def sweep_candidate_windows(window_table, heterozygosity_column="mean_heterozygosity", lower_quantile=0.05):
    """Flag low-diversity windows using a within-table quantile heuristic."""
    if not 0<lower_quantile<1: raise ValueError("lower_quantile must be between 0 and 1")
    df=pd.DataFrame(window_table).copy()
    vals=df[heterozygosity_column].dropna()
    threshold=float(vals.quantile(lower_quantile)) if len(vals) else np.nan
    df["sweep_candidate"]=df[heterozygosity_column]<=threshold
    df["screening_threshold"]=threshold
    return df
