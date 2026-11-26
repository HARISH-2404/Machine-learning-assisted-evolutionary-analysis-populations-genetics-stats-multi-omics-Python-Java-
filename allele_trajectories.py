"""Forward-time trajectory combining selection, mutation and optional drift."""
import numpy as np
from ._core import prob,selected_frequency,ne_int
from .fitness_models import genotype_fitness
from .mutation import mutation_step

def simulate_trajectory(initial_frequency,generations,effective_population_size=None,selection_coefficient=0.0,dominance=0.5,forward_mutation=0.0,reverse_mutation=0.0,fitness=None,seed=None):
    """Order per generation: viability selection -> mutation -> optional WF sampling.

    With Ne=None, trajectory is deterministic. Fitness may be an explicit triplet
    to model overdominance/underdominance; otherwise s,h parameterization is used.
    """
    if int(generations)!=generations or generations<0: raise ValueError('generations must be nonnegative integer')
    p=prob(initial_frequency); f=genotype_fitness(selection_coefficient,dominance) if fitness is None else np.asarray(fitness,float)
    rng=np.random.default_rng(seed); out=[(0,p)]
    ne=None if effective_population_size is None else ne_int(effective_population_size)
    for t in range(1,int(generations)+1):
        p,_=selected_frequency(p,f); p=mutation_step(p,forward_mutation,reverse_mutation)
        if ne is not None: p=rng.binomial(2*ne,p)/(2*ne)
        out.append((t,p))
    return np.asarray(out,float)
