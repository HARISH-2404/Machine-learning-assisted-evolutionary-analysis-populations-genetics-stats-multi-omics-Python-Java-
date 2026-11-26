"""Environment-dependent selection and genotype-by-environment effects."""
import pandas as pd

def environment_selection_schedule(environment_table, generation_col="generation",
                                   environment_col="environment", selection_col="selection_coefficient"):
    required={generation_col,environment_col,selection_col}
    if not required.issubset(environment_table.columns):
        raise ValueError(f"Missing columns: {sorted(required-set(environment_table.columns))}")
    return {(int(r[generation_col]),str(r[environment_col])):float(r[selection_col])
            for _,r in environment_table.iterrows()}

def predict_environmental_trait(genetic_value, environment_effect, gxe_effect=0.0):
    """Simple trait expectation: G + E + GxE; components must be on same scale."""
    return float(genetic_value+environment_effect+gxe_effect)

def apply_environmental_selection(base_s, environment, gxe_coefficients=None):
    """Return s adjusted by an environment-specific effect."""
    return float(base_s+(gxe_coefficients or {}).get(environment,0.0))
