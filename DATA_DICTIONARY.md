# Data dictionary

## Experimental design
`field_trial_design_250x5x3.csv`: plot_id, genotype_id, location_id, replication, block_id. Expected rows: 3,750.

## Genotype dosage
`genotype_dosage_250x300.csv`: genotype_id plus 300 SNP dosage columns. 0 = reference homozygote, 1 = heterozygote, 2 = alternate homozygote; blank = missing. One row per genotype (250 rows).

## Plot phenotype
`plot_phenotypes_3750_plots.csv`: plot_id, genotype_id, location_id, replication, block_id, season and seven phenotype columns. One row per experimental plot (3,750 rows).

## Location environment
`location_environment.csv`: one row per location with illustrative latitude, longitude, altitude, season, temperature, rainfall, soil pH, soil organic carbon and irrigation category.

## Marker map
`marker_map_300.csv`: marker_id, chromosome, physical position, reference/alternate alleles, simulated MAF, observed MAF, missingness and simulated causal labels. Causal labels are truth for simulation evaluation only.

## Genotype metadata and pedigree
- `genotype_information.csv`: genotype ID, genetic group, breeding status and pedigree ID.
- `genotype_sample_metadata.csv`: genotype/sample IDs, species and simplified ploidy note.
- `pedigree_250_genotypes.csv`: genotype and two parent identifiers. Pedigree is illustrative, not a validated pedigree.

## Expression
`expression_250x100.csv`: genotype_id plus 100 continuous synthetic expression features. Not raw counts.

## Processed summaries
- `marker_qc_summary.csv`: marker call rate, missingness, allele frequency, MAF and observed heterozygosity.
- `genotype_dosage_imputed_for_demo.csv`: mean-imputed dosage matrix for demonstration only.
- `genotype_location_overall_means.csv`: genotype means across locations and replications.
- `genotype_by_location_means.csv`: genotype-by-location means across replications.

## Trait definitions
See `reference/annotation/trait_dictionary.csv`. Disease score is simulated on a 1–9 scale, with lower values representing less disease.
