"""Temporal Ne approximation and expected neutral heterozygosity loss."""

def estimate_ne_temporal(p0,p1,generations=1):
    """Approximate Ne=t/(2F), F=(p1-p0)^2/[p0(1-p0)].

    Single-locus estimate is noisy and ignores sampling error, selection and migration.
    """
    p0=float(p0); p1=float(p1); t=int(generations)
    if not 0<p0<1 or not 0<=p1<=1 or t<1: raise ValueError('require 0<p0<1, 0<=p1<=1, generations>=1')
    f=(p1-p0)**2/(p0*(1-p0))
    return float('inf') if f==0 else t/(2*f)

def expected_heterozygosity_decay(H0,effective_population_size,generations):
    """Diploid neutral expectation Ht=H0*(1-1/(2Ne))**t."""
    h=float(H0); ne=int(effective_population_size); t=int(generations)
    if not 0<=h<=1 or ne<1 or t<0: raise ValueError('invalid H0, Ne or generations')
    return h*(1-1/(2*ne))**t
