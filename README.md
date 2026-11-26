# Data package: 250 genotypes × 5 locations × 3 replications × 300 markers

## Experimental design
- Genotypes: 250
- Locations: 5
- Replications per genotype per location: 3
- Total plot records: 3750
- Markers: 300
- Traits: 7
- Molecular expression features: 100

This is a **synthetic demonstration dataset**, not measured experimental data.

## Main files
- `simulated/genotype/genotype_dosage_250x300.csv`: one row per genotype, 300 marker columns.
- `simulated/phenotype/plot_phenotypes_3750_plots.csv`: one row per genotype-location-replication plot.
- `simulated/metadata/field_trial_design_250x5x3.csv`: complete trial layout.
- `simulated/environment/location_environment.csv`: five example location covariate records.
- `simulated/metadata/genotype_sample_metadata.csv`: genotype/sample descriptors.
- `simulated/pedigree/pedigree_250_genotypes.csv`: illustrative pedigree.
- `simulated/expression/expression_250x100.csv`: synthetic continuous expression features.
- `reference/marker_map/marker_map_300.csv`: marker positions and allele labels.
- `reference/annotation/synthetic_gene_annotation.csv`: synthetic gene intervals only.
- `processed/`: marker QC and summarized phenotype tables.

## Experimental design interpretation
Each genotype is evaluated in every location with three replications, giving a balanced multi-environment trial with 3,750 plots. Marker genotype is measured at genotype level and is therefore not repeated as independent DNA observations for each plot.

## Genotype coding and cotton caveat
The 300-marker matrix is simulated as biallelic diploid dosage (0/1/2; missing values are blank). *Gossypium hirsutum* is allotetraploid. For real cotton data, establish whether markers are homeolog-specific, how allele dosage is called, and whether a diploid coding is valid before applying diploid population-genetic estimators.

## Phenotype simulation
Plot phenotypes contain genotype genetic effects, location effects, genotype-by-environment effects, replication/block variation and residual noise. The trait values and units are illustrative; replace them with actual measurement protocols for real experiments.

## Expression data
Expression matrix is synthetic continuous normalized-like data, **not raw RNA-seq counts**. Do not feed it into count-based DESeq2 as if it were raw integer counts.

## Reproducibility
Random seed: 20261127. Files are generated using NumPy/Pandas. `simulated_causal` in marker metadata is simulation truth only and must not be used as a predictor in ML.


## Directory population status

All requested data subdirectories contain data files or README instructions. `raw/` contains synthetic raw-like examples solely to test ingestion; these are not empirical observations. `reference/genome/synthetic_demo_reference.fa` is a toy random sequence for parser testing only, not a cotton reference assembly. Replace these with real source data and a validated *G. hirsutum* reference before biological analysis.
