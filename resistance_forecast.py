"""Candidate resistance-associated allele forecasting and evidence labels."""
import pandas as pd
from .wright_fisher import simulate_wright_fisher, summarize_simulations
from .uncertainty import probability_above, fixation_loss_probability

def forecast_resistance_allele(initial_frequency, generations, population_size,
                               selection_coefficient, dominance=0.5, replicates=1000,
                               threshold=0.8, validated_resistance_evidence=False,
                               seed=1, **kwargs):
    paths=simulate_wright_fisher(initial_frequency,generations,population_size,
        selection_coefficient=selection_coefficient,dominance=dominance,
        replicates=replicates,seed=seed,**kwargs)
    summary=summarize_simulations(paths)
    final=paths[paths.generation==generations]
    return {
      "simulation_paths":paths,
      "frequency_summary":summary,
      "probability_above_threshold":probability_above(paths,threshold,generations),
      "fixation_loss":fixation_loss_probability(paths,generations),
      "evidence_status":"validated_resistance_associated" if validated_resistance_evidence else "candidate_only_not_validated",
      "interpretation":"Predicted allele-frequency dynamics under stated model; not proof of resistance phenotype"
    }
