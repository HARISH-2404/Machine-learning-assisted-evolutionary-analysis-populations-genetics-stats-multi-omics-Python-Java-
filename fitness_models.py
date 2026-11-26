"""Diploid fitness models for biallelic alleles."""
def genotype_fitness(selection_coefficient, dominance=0.5, model="additive"):
    """Return relative fitness (AA, Aa, aa), where A is focal allele.

    Additive: h=0.5; dominant beneficial: h=1; recessive beneficial: h=0.
    Overdominance/underdominance should use explicit genotype_fitnesses.
    """
    s=float(selection_coefficient)
    h={"additive":0.5,"dominant":1.0,"recessive":0.0}.get(model,dominance)
    if model not in {"additive","dominant","recessive","custom"}:
        raise ValueError("model must be additive, dominant, recessive or custom")
    if model=="custom": h=float(dominance)
    w=(1+s,1+h*s,1.0)
    if min(w)<0: raise ValueError("Fitness values cannot be negative")
    return w

def custom_genotype_fitness(w_AA,w_Aa,w_aa):
    """Validate explicit fitnesses; supports overdominance and underdominance."""
    w=tuple(map(float,(w_AA,w_Aa,w_aa)))
    if min(w)<0 or max(w)==0: raise ValueError("Fitnesses must be nonnegative and not all zero")
    return w
