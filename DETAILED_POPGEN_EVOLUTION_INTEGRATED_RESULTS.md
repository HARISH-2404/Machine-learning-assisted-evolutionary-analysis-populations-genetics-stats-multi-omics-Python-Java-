# Detailed Population Genetics, Evolutionary Genetics & Integrated Framework Results

## 1. Executive interpretation

This package analyses the synthetic dataset supplied with the IPGML framework. It contains 250 genotype profiles, 300 biallelic marker columns, 3,750 plot-level observations across 5 locations and three replications per genotype-location combination, and 100 continuous expression features. It includes four genetic groups.

**Key integrated example: SNP_111 on Chr06 (25,226,000 bp).** The alternate allele A is associated with a higher disease score, while the trait dictionary explicitly defines a lower disease score as less disease. Therefore A is labelled *susceptibility-associated* and the reference allele G is labelled *candidate resistance-associated* for this simulation. The association is strong in this synthetic dataset (beta = 1.3044 disease-score units per alternate-allele dosage; R² = 0.3390; p = 8.136e-24; BH q = 2.441e-21). The marker metadata also flags SNP_111 as simulated causal for disease score; this is simulation ground truth, not evidence for a real cotton locus.

The observed allele frequencies are:
- Susceptibility-associated A: **0.0729**
- Candidate resistance-associated G: **0.9271**

This is a real cross-sectional result from the synthetic snapshot plus a separate, clearly marked forward simulation of future generations. The dataset itself does not contain longitudinal allele-count observations.

## 2. Dataset, scope and safeguards

| Feature | Value |
|---|---:|
| Genotypes | 250 |
| Markers | 300 |
| Plot records | 3,750 |
| Locations | 5 |
| Replications per genotype-location | 3 |
| Traits | 7 |
| Expression features | 100 |
| Genetic groups | 4 |
| Random seed used for stochastic simulations | 20261127 |

The genotype matrix is coded as diploid biallelic dosage 0/1/2. *Gossypium hirsutum* is allotetraploid. The current dataset is a diploid demonstration and does not validate cotton homeolog-specific or polyploid dosage estimates. Expression values are continuous synthetic values, not raw RNA-seq counts. The reference sequence and gene coordinates are synthetic demonstration resources.

## 3. Population genetics — results for each listed module

### 3.1 Allele frequency (`allele_frequency.py`)
Calculated reference and alternate allele frequencies, MAF, major/minor allele labels and marker missingness for all 300 markers. Mean MAF = 0.2736; median MAF = 0.2717. Full file: `tables/population_genetics/allele_frequency.csv`. Group-specific frequencies are in `group_allele_frequencies_fis.csv`.

### 3.2 Genotype frequency (`genotype_frequency.py`)
Calculated counts and proportions of 0/0, 0/1 and 1/1 dosage classes for every marker (diploid coding). Full file: `genotype_frequency.csv`.

### 3.3 Hardy–Weinberg equilibrium (`hwe.py`)
Per-marker Pearson chi-square tests were run where all three expected genotype counts were at least five. 239 of 300 markers met that rule. 12 were nominally below p<0.05; **0 remained below BH q<0.05**. SNP_111 itself was not HWE-testable by the expected-count rule because its minor-allele homozygote expectation is too small. See `hwe.csv`; non-testable rows are retained as NA, not treated as evidence of equilibrium.

### 3.4 Heterozygosity (`heterozygosity.py`)
Mean observed heterozygosity Ho = 0.3635; mean expected heterozygosity He = 0.3649. Marker-wise Ho, He and Ho−He are in `heterozygosity.csv`.

### 3.5 Diversity (`diversity.py`)
Calculated biallelic gene diversity (2pq), effective allele number (1/(p²+q²)) and Shannon diversity. Mean gene diversity = 0.3649. Full marker-level values are in `diversity.csv`.

### 3.6 PIC (`pic.py`)
Calculated biallelic polymorphic information content using PIC = 2pq(1−pq). Mean PIC = 0.2913. See `pic.csv`.

### 3.7 FST (`fst.py`)
Pairwise group differentiation was estimated using a Nei-style Ht/Hs ratio from group allele frequencies. The six group-pair mean estimates are in `fst_pairwise_groups.csv`; mean across pairwise group summaries = 0.0155. This is a simple uncorrected proxy for this synthetic dataset, not a full finite-sample corrected Weir–Cockerham estimator.

### 3.8 FIS (`fis.py`)
Calculated per-marker and per-genetic-group FIS = 1−Ho/He. Overall mean marker FIS = 0.0067. Positive values indicate heterozygote deficit; negative values indicate excess under the estimator. See `fis.csv` and `group_allele_frequencies_fis.csv`.

### 3.9 Linkage disequilibrium (`linkage_disequilibrium.py`)
Dosage Pearson correlation and r² were calculated for 44,850 marker pairs. Mean r² = 0.00515; median r² = 0.00233; maximum r² = 0.0942 (SNP_179–SNP_254). This is unphased dosage-correlation LD, not phased D′ or haplotype-based LD. Full output: `linkage_disequilibrium_all_pairs.csv`.

### 3.10 Haplotype (`haplotype.py`)
A six-marker Chr06 haplotype demonstration was generated for SNP_105, SNP_108, SNP_111, SNP_112, SNP_113 and SNP_115. The source genotype data are unphased, so heterozygotes were randomly assigned phase. The resulting haplotype frequencies are **simulated phase illustrations only**, not observed phased haplotypes. See `haplotype_frequencies_simulated_phase.csv` and the method note.

### 3.11 Genetic distance (`genetic_distance.py`)
Calculated pairwise Nei genetic identity/distance between the four genetic groups and standardized Euclidean distances among all genotypes. The latter is a descriptive dosage distance, not a pedigree or kinship estimator. See `genetic_distance_groups.csv` and `genetic_distance_genotypes_euclidean.csv`.

### 3.12 AMOVA (`amova.py`)
A one-factor distance-based AMOVA proxy was calculated using standardized dosage Euclidean distances, grouping samples by `genetic_group`. Estimated PhiPT proxy = 0.0437. This is an exploratory distance-partitioning summary, not a hierarchical AMOVA with permutations or a ploidy-aware genetic distance model. See `amova_group_proxy.csv`.

### 3.13 Population structure (`population_structure.py`)
PCA of standardized mean-imputed marker dosages: PC1 explains 3.96% and PC1+PC2 explain 5.36%. K-means was evaluated from K=2 to 5; highest silhouette was at K=2 (score 0.0316). The low silhouette indicates weak cluster separation. PCA/K-means do not replace STRUCTURE/ADMIXTURE or model-based ancestry estimation. See PCA and cluster result tables.

## 4. Evolutionary genetics — results for each listed module

### 4.1 Selection (`selection.py`)
Cross-sectional genotype-mean association was used to identify allele direction, **not to estimate a selection coefficient**. SNP_111 ALT A dosage is positively associated with disease score. As lower score indicates less disease, A is the susceptibility-associated allele and REF G is candidate resistance-associated. For forward scenarios, assumed selection coefficients s=(0.02, 0.05, 0.1) were tested with dominance h=(0, 0.5, 1). The illustrative primary scenario uses s=0.05 and h=0.5, with relative fitnesses RR=1, RG=1−hs, AA=1−s. See `resistance_candidate_evidence.csv` and forecast tables.

### 4.2 Drift (`drift.py`)
Neutral Wright–Fisher drift was simulated from the observed susceptible ALT frequency across diploid N=25, 50, 100, 250 and 500, with 5,000 replicates and 30 generations. At N=100, the probability that the susceptible ALT was lost by generation 30 was 0.382 in the neutral simulation. See `genetic_drift_neutral_scenarios.csv`.

### 4.3 Mutation (`mutation.py`)
Mutation-only deterministic scenarios were run for susceptible-to-resistance mutation rates μ=(1e-08, 1e-05, 0.0001, 0.001) and reverse rates ν=(0, 1e-08, 1e-05), over 100 generations. These are illustrative parameter settings, not mutation-rate estimates for cotton. See `mutation_scenarios.csv`.

### 4.4 Migration (`migration.py`)
Using SNP_111 ALT frequencies among genetic groups, GROUP_C had the highest frequency (0.1000) and GROUP_A the lowest (0.0403). Migration-only scenarios with m=0, 0.01, 0.05, 0.10 and 0.20 modelled gene flow from the low-susceptibility-allele-frequency group into the high-frequency group. This uses genetic groups as illustrative populations; it does not establish actual migration between them. See `migration_scenarios.csv`.

### 4.5 Recombination (`recombination.py`)
LD decay was simulated as D(t)=D(0)(1−r)^t for recombination fractions r=0, 0.001, 0.01, 0.05, 0.10, 0.25 and 0.50 over 50 generations, initialized from the highest dosage-LD pair. Because the source has no phased haplotypes or genetic map, this is a conceptual decay simulation rather than a recombination-rate estimate. See `recombination_ld_decay_scenarios.csv`.

### 4.6 Selection signatures (`selection_signatures.py`)
A screening score ranked markers using maximum pairwise group FST proxy and low expected heterozygosity. This is an exploratory differentiation/diversity screen. No Tajima's D, SFS-based test, iHS, XP-EHH or demographic null model was inferred because suitable sequence, phased haplotype and longitudinal data are unavailable. See `selection_signature_screen.csv` and top-15 list.

### 4.7 Temporal analysis (`temporal_analysis.py`)
There is only one observed genotype/allele-frequency snapshot (generation not repeated). Therefore observed temporal allele-frequency change, temporal selection coefficient and temporal Ne cannot be estimated. The generation 1–30 rows are labelled model-generated forecasts. `temporal_analysis_simulated_series.csv` records this distinction.

### 4.8 Effective population size (`effective_population_size.py`)
Temporal variance Ne is not estimable from the actual one-timepoint dataset. Instead, a method-recovery simulation used known diploid Ne values 25–1000 and 20,000 replicate one-generation samples. Estimates closely recovered those known simulation values (see `effective_population_size_simulation_recovery.csv`). These are simulation validation outputs, not estimates of the study population's Ne.

## 5. Integrated susceptible-to-resistant allele example

### 5.1 Observed evidence at SNP_111
- Locus: Chr06:25,226,000; REF=G; ALT=A.
- ALT A frequency = 0.0729; REF G frequency = 0.9271.
- Disease-score beta per ALT dosage = 1.3044; R²=0.3390; p=8.136e-24; q=2.441e-21.
- Mean disease score by ALT dosage: dosage 0 = 2.822; dosage 1 = 4.124; dosage 2 = 5.468. Dosage 2 has only one genotype in the dataset.
- Allele frequencies by group are in `resistance_allele_frequency_by_group.csv`.

### 5.2 Generation-by-generation forward simulation
Primary scenario: diploid N=100, s=0.05 against the susceptibility-associated allele, h=0.5, 2,000 Wright–Fisher replicates and 30 generations. It starts from the observed A frequency 0.0729; G starts at 0.9271. At generation 30:
- Mean susceptible ALT A frequency = 0.0390.
- Mean candidate resistance REF G frequency = 0.9610.
- Probability G reached fixation in the simulated replicates = 0.501.
- Probability A was lost = 0.501.

For comparison, neutral drift at N=100 gave probability of A loss = 0.382 by generation 30. The higher loss/fixation probability under selection illustrates the added directional pressure in this chosen scenario; it does **not** prove that the actual cotton population will evolve at this rate. At generation 30 the 95% simulation interval for the susceptible allele frequency is [0.000, 0.210].

### 5.3 Integrated molecular evidence
The integrated file joins marker-trait association, allele frequency, heterozygosity, pairwise group FST screen, nearest synthetic gene annotation and expression associations. SNP_111's nearest annotated synthetic gene is GH_GENE_0038 at 1,990,000 bp. The SNP_111 dosage–nearest-gene expression correlation is -0.010 (p=0.871); nearest-gene expression–disease correlation is 0.032 (p=0.620). There is no expression-based support for this particular locus in this synthetic example. Nearest gene distance is large, and gene annotations are synthetic.

## 6. Files to start with

1. `tables/evolutionary_genetics/resistance_candidate_evidence.csv` — candidate locus evidence and allele labels.
2. `tables/integrated_framework/resistance_allele_frequency_by_group.csv` — current susceptible/resistance-associated frequencies in each group.
3. `tables/forecasting/susceptible_to_resistant_generation_trajectories.csv` — generation-wise forecast means, uncertainty intervals and fixation/loss probabilities across scenarios.
4. `tables/forecasting/primary_susceptible_to_resistant_replicate_paths.csv` — all 2,000 replicate paths for the primary scenario.
5. `tables/integrated_framework/marker_population_genomics_expression_trait_integration.csv` — joined evidence across data layers.
6. `figures/susceptible_to_resistant_trajectory.png` — key allele-transition visualization.

## 7. Scientific limitations and what new data are needed

- For real observed allele changes: genotype/allele counts from the same populations across at least two generations or sampling dates, with sample sizes and denominators.
- For cotton: A/D subgenome-specific marker assignment, verified homeolog specificity, allele dosage and proper tetrasomic/disomic inheritance model.
- For real haplotypes/LD/recombination: phased genotype data and a genetic map or recombination estimates.
- For temporal Ne/migration: repeated sampling and neutral loci, plus sampling design and population connectivity.
- For formal selection signatures: sequence data, suitable demographic nulls, phased haplotypes and appropriate tests.
- For RNA-based integration: raw integer RNA-seq counts, biological replicates, matched genotypes and appropriate gene/homeolog annotations.
- For causal resistance claims: phenotyping under controlled disease challenge, independent population validation, and functional validation.

**Bottom line:** the supplied synthetic data support a strong illustrative association at SNP_111 and permit a transparent model-based scenario in which the resistance-associated REF G frequency increases from 0.927 to a mean of 0.961 over 30 generations under the stated assumptions. The temporal change is simulated, not observed, and must not be reported as a measured evolutionary transition.
