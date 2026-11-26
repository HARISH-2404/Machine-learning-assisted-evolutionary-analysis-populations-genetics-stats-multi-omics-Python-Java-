import numpy as np
from popgenml.evolutionary_genetics import *
assert np.allclose(genotype_fitness(.1,.5),[1.1,1.05,1.0])
assert overdominance_fitness(1,1.2,.9)[1] > 1
assert underdominance_fitness(1,.8,1.1)[1] < 1
p2,mean=selection_step(.2,genotype_fitness(.1,.5)); assert p2>.2 and mean>1
assert simulate_selection(.2,5,.1).shape==(6,3)
d=simulate_drift(.4,10,100,replicates=4,seed=3); assert d.shape==(4,11)
assert 0<=mutation_step(.5,.01,.02)<=1
m=simulate_migration([.1,.9],[[.9,.1],[.2,.8]],3); assert m.shape==(4,2)
t=simulate_trajectory(.2,10,100,selection_coefficient=.05,seed=4); assert t.shape==(11,2)
assert expected_heterozygosity_decay(.5,100,10)<.5
assert estimate_ne_temporal(.4,.41,5)>0
counts=np.array([[8,2],[7,3],[9,1],[6,4]],float)
assert 0<=nucleotide_diversity(counts)<=1
assert np.isfinite(tajimas_d(counts))
print('PASS: all 8 evolutionary genetics modules')
