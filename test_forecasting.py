from popgenml.forecasting import simulate_wright_fisher, summarize_simulations, deterministic_diploid_selection, forecast_allele_and_trait

def test_simulation_shape_and_bounds():
    x=simulate_wright_fisher(.2,5,50,selection_coefficient=.02,replicates=8,seed=7)
    assert len(x)==8*6 and x.allele_frequency.between(0,1).all()

def test_selection_increases_beneficial_allele():
    assert deterministic_diploid_selection(.4,(1.1,1.05,1.0))>.4

def test_forecast_outputs():
    x=forecast_allele_and_trait(.1,3,100,.03,20,1.2,replicates=10,seed=2)
    assert set(["allele_summary","trait_summary","simulation_paths"]).issubset(x)
    assert len(x["allele_summary"])==4
