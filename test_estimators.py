import pandas as pd
from popgenml.integration import estimate_allele_frequency, deterministic_selection_step, estimate_evolutionary_parameters

def sample_counts():
 return pd.DataFrame([
 ["Chr1",100,"G1","A","P1","S1","T1",1,2,10],
 ["Chr1",100,"G1","A","P1","S1","T2",2,5,10],
 ["Chr1",100,"G1","A","P2","S2","T1",1,1,10],
 ["Chr1",100,"G1","A","P2","S2","T2",2,2,10]],
 columns=["chromosome","position","gene_id","allele","population","location","timepoint","generation","allele_count","total_alleles"])

def test_frequency():
 x=estimate_allele_frequency(sample_counts()); assert abs(x.allele_frequency.iloc[0]-.2)<1e-9

def test_selection_update():
 assert deterministic_selection_step(.5,.1)>.5

def test_pipeline_outputs():
 out=estimate_evolutionary_parameters(sample_counts())
 assert "allele_trajectory" in out and "selection_trend" in out
