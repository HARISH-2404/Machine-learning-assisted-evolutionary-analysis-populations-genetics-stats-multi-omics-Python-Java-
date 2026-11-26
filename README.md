# Simulation modules

This package generates reproducible synthetic data for software testing and
method development in population genomics, quantitative genetics and breeding.

## Included modules

- `genotype_simulator.py`: diploid biallelic SNP dosage, population-specific allele
  frequencies, marker map, sample metadata and marker QC.
- `phenotype_simulator.py`: plot-level multi-environment traits with genotype,
  environment, G×E, block, residual, additive and dominance components.
- `expression_simulator.py`: negative-binomial count matrix and known simulated
  differential-expression truth labels.
- `population_simulator.py`: population membership and allele-frequency summaries.
- `evolutionary_simulator.py`: Wright-Fisher drift with selection, dominance,
  migration and mutation.

## Installation

```bash
pip install numpy pandas
```

From repository root, for a `src/` layout:

```bash
export PYTHONPATH="$PWD/src:$PYTHONPATH"
```

## Full experimental-design example

```python
from pathlib import Path
from popgenml.simulation import (
    GenotypeSimulator, PhenotypeSimulator, ExpressionSimulator,
    PopulationSimulator, EvolutionarySimulator
)

out = Path("results/simulation_demo")
genotype, marker_map, sample_meta = GenotypeSimulator(
    n_genotypes=250, n_markers=300, n_chromosomes=15,
    n_populations=5, missing_rate=0.01, seed=20261127
).simulate(out / "genotype")

phenotype = PhenotypeSimulator(
    n_locations=5, replications=3, seed=20261127
).simulate(genotype, marker_map, out / "phenotype")

# Expression samples here are genotypes for demonstration only.
expression, expression_truth = ExpressionSimulator(
    n_genes=100, seed=20261127
).simulate(genotype.genotype_id, output_dir=out / "expression")

population, frequencies = PopulationSimulator(
    n_populations=5, individuals_per_population=50,
    population_names=[f"POP_{i}" for i in range(1,6)]
).simulate(genotype=genotype.assign(population_id=sample_meta.population_id),
           output_dir=out / "population")

evolution = EvolutionarySimulator(
    generations=50, population_size=500, selection_coefficient=0.05,
    dominance_coefficient=0.5, mutation_rate=1e-6,
    migration_rate=0.01, seed=20261127
).simulate(initial_frequencies=marker_map.simulated_source_maf.to_numpy(),
           population_ids=["POP_1"], marker_ids=marker_map.marker_id.tolist(),
           output_dir=out / "evolution")
```

## Model definitions

### Genotype sampling
At marker allele frequency `p`, diploid dosage follows Hardy–Weinberg
probabilities: `P(0)=(1-p)^2`, `P(1)=2p(1-p)`, `P(2)=p^2`.
Population-specific frequencies create simple population differentiation.
Markers are sampled independently; this does not create realistic LD blocks.

### Phenotype model
`Y_ijk = μ + G_i + E_j + (G×E)_ij + B_jk + ε_ijk`.
The genotype term is decomposed into a marker-based additive score and a
heterozygosity-based dominance score. Variance components are simplified
simulation controls, not estimates from a real trial.

### Expression model
Negative-binomial counts have `Var(Y)=μ+αμ²`. Simulated DE labels are retained
in a separate truth file and must never be used as predictor variables in ML.

### Evolutionary model
Genotype fitnesses are `aa=1`, `Aa=1+h·s`, `AA=1+s`. Here `h=0` represents a
recessive beneficial allele, `h=0.5` an additive fitness effect and `h=1` a
dominant beneficial allele. This model does not include heterozygote advantage
(overdominance, where heterozygote fitness exceeds both homozygotes); that
requires a different three-genotype fitness specification.

## Scientific cautions

All outputs are synthetic. They are not observed field, sequencing or
population data. The genotype model assumes diploid 0/1/2 dosage. Cotton
(*Gossypium hirsutum*) is allotetraploid; real cotton analyses require a
validated ploidy/dosage model and appropriate A/D subgenome and homeolog
handling. The expression simulator makes synthetic counts, not experimental
RNA-seq replicates. The evolutionary model is a simplified independent-locus
Wright-Fisher model and does not simulate haplotypes, recombination, linkage
disequilibrium, epistasis, or polyploid inheritance.
