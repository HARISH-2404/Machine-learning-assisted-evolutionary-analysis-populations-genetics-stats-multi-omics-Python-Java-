# Multi-omics integration

Python package for sample harmonization, matched-block merging, training-only scaling, missingness filtering/imputation, cross-omics correlation screening, and early/late integration.

## Data convention
Each input block is a `pandas.DataFrame`: rows = biological samples, columns = features. Use stable sample IDs and feature IDs. Prefixing feature columns avoids collisions across omics.

## Example
```python
import pandas as pd
from popgenml.integration import merge_omics_blocks, scale_feature_block
blocks = {"geno": genotype_df, "expr": expression_df, "pheno": phenotype_df}
X = merge_omics_blocks(blocks, join="inner")
X_train, X_test = scale_feature_block(X_train, X_test, method="zscore")
```

## Critical safeguards
- Confirm that IDs refer to the same biological individual, tissue, treatment, developmental stage and time point. Do not equate plant, plot, RNA library and genotype IDs without a documented mapping.
- Inspect ID-normalization collisions; the helper raises an error rather than silently merging them.
- Fit imputation, scaling, filtering and feature selection inside each training fold. Never fit using the full dataset before validation.
- Correct batch effects within an appropriate experimental design; do not remove biological signal confounded with batch.
- Omics blocks have different distributions and measurement scales. Normalize each block appropriately before integration; do not compare raw genotype dosage, RNA counts and phenotype values directly.
- Correlation is association, not regulation or causation. Correct for multiple testing in downstream inference.
- For cotton, account for allotetraploid A- and D-subgenomes, homeologous loci, allele dosage/ploidy, and ambiguous read/variant assignment. This toolkit does not perform polyploid-aware genotype calling or homeolog resolution.

## Scope
This is a transparent baseline toolkit, not a full implementation of MOFA+, mixOmics/DIABLO, iCluster, or a causal multi-omics model. Use validated specialist tools for those analyses.


## Added interaction modules

- `genotype_expression.py`: aligns sample-level genotype and expression blocks, optionally aggregates variant features to gene IDs.
- `genotype_phenotype.py`: joins marker data to phenotype observations and optional metadata.
- `multiomics.py`: merges multiple blocks, reports coverage and screens cross-block correlations.
- `evolutionary_integration.py`: tracks an allele at chromosome/position/gene across population, location, time point and generation; calculates allele frequency, frequency change and change per generation; optionally joins quantitative-genetic and phenotype columns.

### Longitudinal input schema

`allele_counts` must contain: `chromosome`, `position`, `gene_id`, `allele`, `population`, `location`, `timepoint`, `generation`, `allele_count`, `total_alleles`.

Optional `trait_data` should include the matching keys `population`, `location`, `timepoint`, `generation` plus traits/parameters such as fiber strength, disease score, yield, boll weight, heritability estimates or selection differential. Include a clear `sample_id` or experimental-unit key in your own analysis if multiple records share the same group keys; avoid accidental many-to-many joins.

The trajectory output includes allele frequency and temporal change. A shift from a susceptible-associated allele to a resistance-associated allele can be displayed only when those allele effects have been independently established (e.g. validated phenotype association/functional evidence). Frequency increase by itself does not establish resistance, selection, or causality. Include untreated/control populations and replicate populations where possible to distinguish selection from drift, migration, bottlenecks and sampling noise.

For diploid biallelic data, `total_alleles` is typically twice the number of successfully genotyped individuals. For allotetraploid cotton, allele dosage and denominator must reflect the specific locus/homeolog and ploidy-aware calling; do not blindly use 2N. Position coordinates must state genome assembly and coordinate convention (0-based vs 1-based).

## Detailed estimation modules

- `allele_frequency_estimation.py`: allele frequency, expected heterozygosity, FIS and generation-wise changes.
- `selection_estimation.py`: approximate logit-frequency slope and deterministic diploid viability-selection simulation with dominance.
- `drift_migration_estimation.py`: rough temporal-Ne, one-source admixture estimate and neutral heterozygosity-decay expectation.
- `population_differentiation.py`: expected heterozygosity, pairwise approximate FST and diversity summaries.
- `quantitative_parameter_estimation.py`: simplified broad-sense heritability, genetic correlation, selection differential, realized response and additive/dominance effects.
- `locus_trait_association.py`: single-marker dosage regression screen with BH-style FDR.
- `evolutionary_estimation_pipeline.py`: orchestrates estimates and joins trait/expression observations.

### Minimum allele-count input
One row per chromosome-position-gene-allele-population-location-time-generation:
`chromosome, position, gene_id, allele, population, location, timepoint, generation, allele_count, total_alleles`.

### Interpretation limits
Selection logit slope is a descriptive approximation, not a universal estimate of diploid selection coefficient. Temporal Ne needs many approximately neutral, independent loci and adequate sampling. FST here is an uncorrected heterozygosity approximation. Heritability is a simplified balanced-design estimator; use REML/mixed models for unbalanced trials, GxE and spatial effects. Marker-trait regression is not mixed-model GWAS. Adjust for structure, kinship, batch and multiple testing. Avoid using the same data to select features and evaluate predictive accuracy.

For cotton, estimate allele counts at the correct A/D homeolog and ploidy; standard diploid 0/1/2 dosage and 2N allele denominator are not universally valid.
