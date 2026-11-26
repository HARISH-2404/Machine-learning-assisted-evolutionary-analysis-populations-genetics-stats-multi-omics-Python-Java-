"""Pearson phenotypic/genetic correlation from paired vectors."""
import math

def _corr(x,y):
    if len(x)!=len(y) or len(x)<3: raise ValueError("Equal vectors with >=3 observations required")
    mx=sum(x)/len(x); my=sum(y)/len(y)
    sx=sum((a-mx)**2 for a in x); sy=sum((b-my)**2 for b in y)
    if sx==0 or sy==0: return float("nan")
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(sx*sy)

def genetic_correlation(trait_x, trait_y):
    """Correlation of genotype means; interpret as genetic only if design supports it."""
    return _corr(list(map(float,trait_x)),list(map(float,trait_y)))
