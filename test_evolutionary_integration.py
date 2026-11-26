import pandas as pd
from popgenml.integration.evolutionary_integration import build_evolutionary_trajectory, summarize_locus_population

def test_frequency_change_over_generations():
    d=pd.DataFrame([
      ["Chr1",100,"G1","A", "P1","Site1","T1",1,2,10],
      ["Chr1",100,"G1","A", "P1","Site1","T2",2,8,10]],
      columns=["chromosome","position","gene_id","allele","population","location","timepoint","generation","allele_count","total_alleles"])
    out=build_evolutionary_trajectory(d)
    assert out.allele_frequency.tolist()==[0.2,0.8]
    assert abs(out.frequency_change.iloc[1]-0.6)<1e-9
    assert abs(summarize_locus_population(out).net_frequency_change.iloc[0]-0.6)<1e-9
