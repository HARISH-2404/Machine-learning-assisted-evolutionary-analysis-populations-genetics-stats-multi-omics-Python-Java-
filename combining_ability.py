"""Diallel combining-ability summaries from a cross-mean table."""
def gca_sca(cross_means):
    """Balanced diallel approximation.

    cross_means: {(parent_i,parent_j): mean}. Estimates GCA as row/column
    marginal deviations from grand mean; SCA is cross deviation after GCA.
    Reciprocal effects and incomplete/unbalanced designs are not modelled.
    """
    if not cross_means: raise ValueError("cross_means cannot be empty")
    parents=sorted({p for pair in cross_means for p in pair})
    grand=sum(map(float,cross_means.values()))/len(cross_means)
    effects={}
    for p in parents:
        vals=[float(v) for (a,b),v in cross_means.items() if a==p or b==p]
        if not vals: effects[p]=0.0
        else: effects[p]=sum(vals)/len(vals)-grand
    sca={}
    for (a,b),v in cross_means.items():
        sca[(a,b)]=float(v)-grand-effects[a]-effects[b]
    return {"grand_mean":grand,"gca":effects,"sca":sca}
