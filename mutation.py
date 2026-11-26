"""Two-way mutation at a biallelic locus."""
from ._core import prob

def mutation_step(p,forward_rate,reverse_rate=0.0):
    """u is A->a and v is a->A: p' = p(1-u)+(1-p)v."""
    p=prob(p,'p'); u=prob(forward_rate,'forward_rate'); v=prob(reverse_rate,'reverse_rate')
    return p*(1-u)+(1-p)*v

def simulate_mutation(initial_frequency,generations,forward_rate,reverse_rate=0.0):
    if int(generations)!=generations or generations<0: raise ValueError('generations must be nonnegative integer')
    p=prob(initial_frequency); out=[(0,p)]
    for t in range(1,int(generations)+1): p=mutation_step(p,forward_rate,reverse_rate); out.append((t,p))
    return out
