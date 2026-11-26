"""Pairwise dosage LD (Pearson r and r²); optional chromosome/window filters.

This estimates dosage correlation, not phased-haplotype D-prime. For each
pair, only individuals called at both loci are used.
"""
import numpy as np
import pandas as pd
from itertools import combinations
from ._core import validate_dosage,save

def linkage_disequilibrium(genotypes,marker_map=None,id_col="genotype_id",
                           min_complete=20,max_distance_bp=None,output=None):
    markers=validate_dosage(genotypes,id_col)
    if min_complete<3: raise ValueError("min_complete must be >=3.")
    mp={}
    if marker_map is not None:
        if not {"marker_id","position_bp"}.issubset(marker_map): raise ValueError("marker_map requires marker_id and position_bp.")
        mp=marker_map.set_index("marker_id").to_dict("index")
    x=genotypes[markers].to_numpy(float); rows=[]
    for i,j in combinations(range(len(markers)),2):
        a,b=markers[i],markers[j]
        d=np.nan
        if a in mp and b in mp and mp[a].get("chromosome")==mp[b].get("chromosome"):
            d=abs(float(mp[a]["position_bp"])-float(mp[b]["position_bp"]))
        if max_distance_bp is not None and np.isfinite(d) and d>max_distance_bp: continue
        valid=np.isfinite(x[:,i])&np.isfinite(x[:,j]); n=int(valid.sum())
        if n<min_complete: continue
        u,v=x[valid,i],x[valid,j]
        if np.std(u)==0 or np.std(v)==0: continue
        r=float(np.corrcoef(u,v)[0,1])
        rows.append({"marker_1":a,"marker_2":b,"n_complete":n,
                     "distance_bp":d,"r":r,"r_squared":r*r})
    out=pd.DataFrame(rows,columns=["marker_1","marker_2","n_complete","distance_bp","r","r_squared"])
    if output: save(out,output)
    return out
