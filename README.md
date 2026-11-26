# Configuration guide

These YAML files centralize analysis parameters. Python modules should load and validate them rather than hard-code thresholds.

- `simulation.yaml`: population, genotype, phenotype, expression and evolutionary simulation.
- `population_genetics.yaml`: QC, allele/genotype frequencies, HWE, diversity, PIC, F-statistics, LD, haplotypes, distance, structure and AMOVA.
- `genomics.yaml`: variant filtering, GWAS, QTL, annotation, selection scans and genomic relationship matrices.
- `transcriptomics.yaml`: RNA-seq QC, normalization, differential expression, eQTL, TWAS, co-expression and allele-specific expression.
- `breeding.yaml`: variance components, heritability, genetic advance, BLUP-family methods, multi-trait analysis and genomic selection.
- `machine_learning.yaml`: preprocessing, grouped validation, model settings, feature selection, explainability and forecasting.

## Implementation contract
A YAML option is a requested configuration, not proof that the corresponding algorithm is implemented. Each executable module must validate the keys it consumes, report unsupported options and record the configuration used. Do not silently ignore enabled methods.

## Reproducibility
Record software versions, input checksums, genome assembly, ploidy assumptions, marker filters and sample exclusions in every analysis report.

## Biological cautions
The generic dosage settings assume biallelic diploid markers coded 0/1/2. Cotton is allotetraploid; homeologous loci, allele dosage and subgenome-specific calling must be handled explicitly before applying diploid assumptions to cotton data.
