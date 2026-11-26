"""Expected genetic advance under truncation selection."""
import statistics

def genetic_advance(mean, phenotypic_sd, heritability, selection_intensity):
    """GA = i * H2 * sigmaP; returns GA and percent of mean."""
    if phenotypic_sd<0 or not 0<=heritability<=1 or selection_intensity<0:
        raise ValueError("Invalid SD, heritability, or selection intensity")
    ga=selection_intensity*heritability*phenotypic_sd
    return {"genetic_advance":ga,"genetic_advance_percent":100*ga/mean if mean else float("nan"),
            "predicted_selected_mean":mean+ga}
