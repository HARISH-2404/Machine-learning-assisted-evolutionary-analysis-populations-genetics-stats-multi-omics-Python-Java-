"""Deterministic viability selection."""
import numpy as np
from ._core import prob, genotype_frequencies, selected_frequency
from .fitness_models import genotype_fitness

def selection_step(p, fitness):
    """Apply selection and return (post-selection allele-A frequency, mean fitness)."""
    return selected_frequency(prob(p,'p'),fitness)

def simulate_selection(initial_frequency,generations,selection_coefficient,dominance=0.5,fitness=None):
    """Return ndarray columns generation,p_A,mean_fitness including generation zero."""
    if int(generations)!=generations or generations<0: raise ValueError('generations must be a nonnegative integer')
    f=genotype_fitness(selection_coefficient,dominance) if fitness is None else np.asarray(fitness,float)
    p=prob(initial_frequency); out=[]
    for t in range(int(generations)+1):
        out.append((t,p,float(genotype_frequencies(p)@f)))
        if t<int(generations): p,_=selection_step(p,f)
    return np.asarray(out,float)
