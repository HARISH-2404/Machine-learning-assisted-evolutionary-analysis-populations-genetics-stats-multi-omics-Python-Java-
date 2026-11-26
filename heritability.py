"""Broad- and narrow-sense heritability calculations."""
def broad_sense_heritability(Vg, Vp):
    if Vp<=0 or Vg<0: raise ValueError("Require Vp > 0 and Vg >= 0")
    return Vg/Vp

def narrow_sense_heritability(Va, Vp):
    if Vp<=0 or Va<0: raise ValueError("Require Vp > 0 and Va >= 0")
    return Va/Vp

def entry_mean_heritability(Vg, Vge, Ve, environments, replications):
    """H2 for genotype means under a balanced multi-environment design."""
    if min(Vg,Vge,Ve)<0 or environments<1 or replications<1:
        raise ValueError("Variances must be nonnegative and counts positive")
    denom=Vg+Vge/environments+Ve/(environments*replications)
    return Vg/denom if denom else float("nan")
