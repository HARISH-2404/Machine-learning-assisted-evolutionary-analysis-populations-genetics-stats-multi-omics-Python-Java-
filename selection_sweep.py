"""Site-frequency-spectrum summaries used in selection-sweep screening."""
import math
import numpy as np

def nucleotide_diversity(allele_counts):
    """Mean per-site pairwise nucleotide diversity from n_sites x 2 allele counts."""
    c=np.asarray(allele_counts,float)
    if c.ndim!=2 or c.shape[1]!=2 or np.any(c<0) or not np.all(np.isfinite(c)): raise ValueError('counts must be finite n_sites x 2 nonnegative values')
    n=c.sum(axis=1); pi=np.zeros(len(n)); ok=n>1
    pi[ok]=2*c[ok,0]*c[ok,1]/(n[ok]*(n[ok]-1))
    return float(pi.mean()) if len(pi) else float('nan')

def tajimas_d(allele_counts):
    """Tajima's D for equal-sized haploid chromosome samples across biallelic sites.

    Input counts each row [ancestral/reference count, alternate count]. Requires
    >=2 segregating sites; returns NaN if variance denominator is zero.
    """
    c=np.asarray(allele_counts,float)
    if c.ndim!=2 or c.shape[1]!=2 or len(c)==0 or np.any(c<0): raise ValueError('counts must be nonempty n_sites x 2')
    n=c.sum(axis=1)
    if np.any(n!=n[0]) or n[0]<2: raise ValueError('all sites must have equal sample size >=2')
    N=int(n[0]); S=int(np.sum((c[:,0]>0)&(c[:,1]>0)))
    if S<2: return float('nan')
    a1=sum(1/i for i in range(1,N)); a2=sum(1/(i*i) for i in range(1,N))
    b1=(N+1)/(3*(N-1)); b2=2*(N*N+N+3)/(9*N*(N-1))
    c1=b1-1/a1; c2=b2-(N+2)/(a1*N)+a2/(a1*a1)
    e1=c1/a1; e2=c2/(a1*a1+a2); theta=S/a1
    pi=float(np.sum(2*c[:,0]*c[:,1]/(N*(N-1))))
    den=math.sqrt(e1*S+e2*S*(S-1))
    return (pi-theta)/den if den>0 else float('nan')
