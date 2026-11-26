# Detailed Results Report — Integrated Population Genomics and Machine Learning

## Dataset and analysis scope
All calculations in this report use the synthetic dataset `250 genotypes × 5 locations × 3 replications × 300 markers`, seed 20261127. There are 3,750 plot observations, 250 genotype-level marker profiles, 300 diploid-coded biallelic markers, 100 synthetic continuous expression features and 7 traits. All outputs are simulation results, not empirical cotton findings.

## 1. Data structure and quality
- Genotype matrix: 250 × 300.
- Phenotype matrix: 3750 plots, balanced 5 locations and 3 replications per genotype-location cell.
- Mean genotype-marker missingness: 1.0293%.
- Expression matrix: 250 genotypes × 100 features; values are continuous synthetic normalized-like values, not raw RNA-seq counts.
- The genotype matrix uses diploid 0/1/2 dosage. *Gossypium hirsutum* is allotetraploid, so these data do not validate a polyploid cotton analysis.

## 2. Population genetics
Per-marker outputs are in `tables/population_genetics/marker_population_genetics.csv`. It includes called sample count, missingness, alternate allele frequency, MAF, observed/expected heterozygosity, biallelic gene diversity/PIC, HWE chi-square and p-value, chromosome and position.
- Mean MAF: 0.2736; median MAF: 0.2717.
- Mean observed heterozygosity: 0.3635.
- Mean expected heterozygosity: 0.3649.
- Genotype frequencies are in `genotype_frequencies.csv`.
- All pairwise dosage correlations (r and r²) are in `linkage_disequilibrium_pairs.csv`; 44,850 marker pairs were evaluated, with mean r² 0.0053. This is dosage-correlation LD, not phased haplotype D′.
- Group differentiation proxy by genetic group is in `group_differentiation_proxy.csv`. Treat it as exploratory; it is not a fully corrected population FST estimator.

## 3. Quantitative genetics and multi-environment traits
Plot-level descriptive statistics are in `trait_descriptive_statistics.csv`; location summaries are in `trait_by_location.csv`; genotype means are in `genotype_trait_means.csv`.
- `variance_components_approx.csv` gives balanced two-way ANOVA method-of-moments estimates for genotype, location, G×E and residual components, plus broad-sense entry-mean H². These are approximations under a balanced random-effects interpretation, not REML estimates.
- `genetic_mean_trait_correlations.csv` contains Pearson correlations among genotype means.
- `selection_differential_top10.csv` reports the difference between the top 10% selected genotype mean and the overall mean. This is a descriptive selection differential, not a predicted response to selection.

## 4. Genomics and marker-trait association
`single_marker_association_all_traits.csv` contains 300-marker linear regressions for each of seven traits (2,100 marker-trait tests), including per-dosage beta, R², p-value and BH q-value. `top20_markers_per_trait.csv` provides ranked marker tables. These are exploratory single-marker scans, not mixed-model GWAS or validated QTL. Population structure and relatedness correction are not applied.

## 5. Transcriptomic / multi-omics integration
`expression_trait_correlations.csv` reports 100 expression-feature × 7 trait correlations with BH adjustment. `top1000_marker_expression_associations.csv` reports the top 1,000 marker-expression correlations. Expression values are synthetic continuous features, not RNA-seq counts; these results do not represent differential expression, eQTL causal effects or regulatory proof.

## 6. Machine learning
`model_comparison_5fold_cv.csv` compares Ridge, Random Forest, Extra Trees, Gradient Boosting and RBF-SVR for each trait using 5-fold shuffled genotype-level cross-validation. The prediction target is the genotype mean phenotype and the features are 300 marker dosages. The best model per trait by mean CV R² is:
| trait                     | model            |   mean_R2 |   mean_RMSE |   mean_MAE |
|:--------------------------|:-----------------|----------:|------------:|-----------:|
| boll_weight_g             | GradientBoosting |  0.737195 |    0.269258 |   0.2066   |
| disease_score             | GradientBoosting |  0.770153 |    0.380911 |   0.295707 |
| fiber_length_mm           | GradientBoosting |  0.769769 |    0.548815 |   0.416545 |
| fiber_strength_g_tex      | GradientBoosting |  0.658163 |    0.952737 |   0.733107 |
| lint_percent              | GradientBoosting |  0.688355 |    0.985648 |   0.774843 |
| micronaire                | GradientBoosting |  0.646389 |    0.220145 |   0.174981 |
| seed_cotton_yield_kg_plot | GradientBoosting |  0.678794 |    0.241546 |   0.193214 |
`top30_random_forest_marker_importance.csv` contains feature-importance rankings from models fit to all genotypes. Feature importance is predictive attribution, not causal evidence. This validation is a baseline; family-aware, population-aware and external validation are still needed.

## 7. Evolutionary genetics and forecasting
The dataset has no repeated allele counts by generation/timepoint. Therefore observed selection coefficients, allele-frequency trajectories through time, temporal effective population size, drift estimates from time series, fixation/loss observed across generations and resistance transition rates cannot be estimated from this dataset.
To demonstrate forecasting code, `hypothetical_wright_fisher_scenarios.csv` contains scenario simulations for the 10 strongest nominal fiber-strength marker associations, starting from observed alternate-allele frequency and applying hypothetical selection coefficients s = 0, 0.02, 0.05, 0.10, dominance h = 0.5, N = 100 diploid individuals, 500 replicates and 20 generations. These are what-if projections only, not inferred biological forecasts. The explicit limitation is also recorded in `temporal_data_limitation.csv`.

## 8. Parameters not estimated from available inputs
The current synthetic data do not contain explicit time-series allele counts, raw RNA-seq counts, phased haplotypes, full sequence-level variant calls, replicated disease challenge time series, or a validated cotton polyploid dosage model. Consequently, robust Tajima's D/SFS-based demography, temporal Ne, migration rate, selection sweep inference, haplotype-based LD, formal DESeq2/edgeR differential expression, ASE, TWAS, causal regulatory networks, polyploid GWAS/QTL, and validated resistance evolution are not reported as empirical estimates.

## Reproducibility and interpretation
Random seed for scenario projections: 20261127. All tables are CSV. Figures are PNG. Use the dataset's `simulated_causal` truth columns only for post hoc validation of simulation recovery, never as predictor variables. No claim in this report should be interpreted as a result from actual field or molecular experiments.
