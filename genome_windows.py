"""Window-based genomic summaries and marker density."""
from collections import defaultdict
import numpy as np

def window_coordinates(length, window_size, step=None):
    """Generate 1-based inclusive sliding windows."""
    length=int(length); window_size=int(window_size); step=window_size if step is None else int(step)
    if length<1 or window_size<1 or step<1: raise ValueError("length, window_size and step must be positive")
    out=[]; start=1
    while start<=length:
        end=min(start+window_size-1,length)
        out.append((start,end))
        if end==length: break
        start+=step
    return out

def summarize_variant_windows(variants, window_size=100000, step=None):
    """Count variants per chromosome/window; variants require chrom and 1-based pos."""
    groups=defaultdict(list)
    for v in variants: groups[str(v["chrom"])].append(int(v["pos"]))
    rows=[]
    for chrom,positions in groups.items():
        maxpos=max(positions)
        for start,end in window_coordinates(maxpos,window_size,step):
            count=sum(start<=p<=end for p in positions)
            rows.append({"chrom":chrom,"start":start,"end":end,"variant_count":count,
                         "density_per_kb":count/((end-start+1)/1000)})
    return rows

def nucleotide_diversity_windows(positions, pi_values, chromosome_length, window_size=100000):
    """Mean per-site pi per window; positions are 1-based."""
    if len(positions)!=len(pi_values): raise ValueError("positions and pi_values lengths differ")
    p=np.asarray(positions,dtype=int); pi=np.asarray(pi_values,dtype=float)
    rows=[]
    for start,end in window_coordinates(chromosome_length,window_size):
        m=(p>=start)&(p<=end)&np.isfinite(pi)
        rows.append({"start":start,"end":end,"n_sites":int(m.sum()),
                     "mean_pi":float(pi[m].mean()) if m.any() else float("nan")})
    return rows
