"""End-to-end forecast orchestration for one focal allele."""
from .wright_fisher import simulate_wright_fisher, summarize_simulations
from .uncertainty import probability_above, fixation_loss_probability
from .trait_forecasting import project_trait_from_allele_frequency

def forecast_allele_and_trait(p0, generations, population_size, selection_coefficient,
                              baseline_trait_mean, additive_effect, dominance=0.5,
                              dominance_effect=0.0, replicates=500, seed=1, **kwargs):
    paths=simulate_wright_fisher(p0,generations,population_size,
        selection_coefficient=selection_coefficient,dominance=dominance,
        replicates=replicates,seed=seed,**kwargs)
    summary=summarize_simulations(paths)
    trait_paths=project_trait_from_allele_frequency(paths,baseline_trait_mean,
        additive_effect,baseline_frequency=p0,dominance_effect=dominance_effect)
    trait_summary=trait_paths.groupby("generation").projected_trait_mean.agg(
        mean="mean",sd="std",median="median").reset_index()
    return {"simulation_paths":paths,"allele_summary":summary,
            "trait_paths":trait_paths,"trait_summary":trait_summary,
            "probability_above_0_8":probability_above(paths,0.8),
            "fixation_loss":fixation_loss_probability(paths)}
