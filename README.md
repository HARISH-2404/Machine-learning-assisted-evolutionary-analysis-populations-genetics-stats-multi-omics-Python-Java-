# Evolutionary forecasting

This package extends population-genetic estimation into forward prediction of allele-frequency trajectories and associated trait means.

## Modules
- `fitness_models.py`: additive, dominant, recessive and custom genotype fitness; custom triplets support overdominance/underdominance.
- `wright_fisher.py`: diploid Wright–Fisher forward simulation with selection, drift, mutation and migration; replicate paths and quantile summaries.
- `environment_gxe.py`: environment-specific selection schedules and simple G×E trait expectation.
- `uncertainty.py`: probability of crossing frequency thresholds, fixation/loss probabilities and bootstrap intervals.
- `trait_forecasting.py`: single-locus additive/dominance and independent-locus polygenic trait projections.
- `resistance_transition.py`: evidence-aware candidate resistance/susceptibility allele trajectory summaries.
- `forecast_pipeline.py`: combined allele and trait forecasting workflow.

## Quick start
```python
from popgenml.forecasting import forecast_allele_and_trait
result = forecast_allele_and_trait(
    p0=0.15, generations=20, population_size=200,
    selection_coefficient=0.05, baseline_trait_mean=25.0,
    additive_effect=1.2, dominance=0.5,
    replicates=1000, seed=2026
)
result["allele_summary"].to_csv("allele_forecast.csv", index=False)
result["trait_summary"].to_csv("trait_forecast.csv", index=False)
```

## Model assumptions and warnings
- The default Wright–Fisher model is a single diploid, biallelic locus with random mating, discrete non-overlapping generations and binomial drift sampling of 2N gene copies.
- Selection fitness is genotype-based: AA=1+s, Aa=1+h*s, aa=1. Custom fitness triplets allow overdominance (wAa greater than both homozygotes) or underdominance.
- Mutation and migration are one-generation frequency transformations. Migration requires an explicitly supplied migrant allele frequency.
- Environment effects in `environment_gxe.py` are helper functions; the caller must supply a generation-specific selection schedule and run environment-specific simulations. The core simulator accepts generation-to-selection adjustments.
- Trait projection assumes Hardy–Weinberg proportions and supplied effects. Single-locus and independent-locus projections are not equivalent to GBLUP or a full quantitative-genetic mixed model.
- A frequency increase does not prove selection or resistance. Validate resistance effects with phenotyping, independent association/functional evidence and control populations.
- Forecast uncertainty reflects stochastic replicates under the specified model, not uncertainty in parameter estimates unless parameters are also sampled.
- Cotton is allotetraploid. Do not apply diploid 2N sampling to cotton loci without homeolog-specific ploidy-aware allele dosage and an appropriate inheritance model. This implementation is a diploid baseline, not a validated cotton polyploid simulator.


## Additional requested forecasting modules

- `allele_trajectory.py`: builds observed locus-specific frequency trajectories from chromosome/position/gene/allele/population/location/generation allele counts and forecasts a focal locus.
- `phenotype_trajectory.py`: summarizes observed trait means and projects a single-locus trait trajectory.
- `resistance_forecast.py`: predicts candidate resistance-associated allele frequency and threshold/fixation probabilities; labels evidence status explicitly.
- `breeding_response.py`: breeder's equation response (R=h²S), realized response and linear cycle projection.
- `future_population.py`: forecasts independent populations with population-specific parameters and converts frequencies to expected allele-copy counts under diploid assumptions.

### Important assumptions
Breeder's equation uses **narrow-sense** heritability (h²), not broad-sense H². Trait projections require supplied effect estimates and assumptions about genotype frequencies. `future_population.py` uses 2N for allele-copy counts and is diploid only. These modules do not yet model structured migration between multiple populations jointly; each population is simulated independently unless the caller explicitly supplies a migration source frequency in the single-population simulator. For cotton, replace diploid dosage and 2N assumptions with validated homeolog-specific polyploid models.
