"""Future population composition and multi-population projections."""
import pandas as pd
from .wright_fisher import simulate_wright_fisher, summarize_simulations

def forecast_multiple_populations(population_parameters, generations, replicates=500, seed=1):
    """Run independent forecasts by population.

    population_parameters maps population label to p0, population_size,
    selection_coefficient and optional mutation/migration/fitness parameters.
    """
    paths=[]; summaries=[]
    for i,(name,params) in enumerate(population_parameters.items()):
        p=dict(params)
        p0=p.pop("p0"); n=p.pop("population_size"); s=p.pop("selection_coefficient")
        sim=simulate_wright_fisher(p0,generations,n,s,replicates=replicates,
                                   seed=seed+i,**p)
        sim.insert(0,"population",name); paths.append(sim)
        sm=summarize_simulations(sim.drop(columns="population"))
        sm.insert(0,"population",name); summaries.append(sm)
    return pd.concat(paths,ignore_index=True),pd.concat(summaries,ignore_index=True)

def project_population_counts(frequency_summary, population_sizes, frequency_col="mean"):
    """Convert predicted frequency to expected focal-allele copy count (diploid 2N)."""
    x=frequency_summary.copy()
    x["population_size"]=x["population"].map(population_sizes)
    if x["population_size"].isna().any(): raise ValueError("Missing population size for one or more populations")
    x["expected_allele_copies"]=2*x["population_size"]*x[frequency_col]
    return x
