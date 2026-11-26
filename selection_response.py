"""Breeder's equation and correlated response helpers."""
def response_to_selection(selection_differential, heritability):
    if not 0<=heritability<=1: raise ValueError("Heritability must be in [0,1]")
    return {"response":heritability*selection_differential,
            "response_per_generation":heritability*selection_differential}

def correlated_response(selection_intensity, genetic_correlation,
                        target_genetic_sd, selected_trait_heritability):
    if not -1<=genetic_correlation<=1 or target_genetic_sd<0 or not 0<=selected_trait_heritability<=1:
        raise ValueError("Invalid correlation, SD or heritability")
    return selection_intensity*genetic_correlation*target_genetic_sd*selected_trait_heritability
