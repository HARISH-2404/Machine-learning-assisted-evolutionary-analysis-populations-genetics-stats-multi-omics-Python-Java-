package org.ipgml.populationgenetics;
import java.nio.file.*;
/** CLI example: java ... PopulationGeneticsCli genotypes.csv output_directory */
public class PopulationGeneticsCli {
 public static void main(String[] args)throws Exception{
  if(args.length<2){System.err.println("Usage: PopulationGeneticsCli genotypes.csv output_dir");System.exit(2);}
  GenotypeData d=GenotypeData.readCsv(Paths.get(args[0]));Path out=Paths.get(args[1]);Files.createDirectories(out);
  GenotypeData.writeCsv(out.resolve("allele_frequency.csv"),
   "population_id,marker_id,n_called,call_rate,alt_frequency,ref_frequency,maf",AlleleFrequency.estimate(d));
  GenotypeData.writeCsv(out.resolve("genotype_frequency.csv"),
   "population_id,marker_id,n_called,n_0,n_1,n_2,freq_0,freq_1,freq_2",GenotypeFrequency.estimate(d));
  GenotypeData.writeCsv(out.resolve("hwe.csv"),
   "population_id,marker_id,n,p_alt,obs_00,obs_01,obs_11,exp_00,exp_01,exp_11,chi_square,p_value,expected_counts_all_ge_min",Hwe.test(d));
  GenotypeData.writeCsv(out.resolve("heterozygosity.csv"),
   "population_id,marker_id,n,Ho,He,He_unbiased",Heterozygosity.estimate(d));
  GenotypeData.writeCsv(out.resolve("diversity.csv"),
   "population_id,marker_id,Ho,He,effective_allele_number",Diversity.estimate(d));
  GenotypeData.writeCsv(out.resolve("pic.csv"),
   "population_id,marker_id,n,p_alt,q_ref,PIC",Pic.estimate(d));
  GenotypeData.writeCsv(out.resolve("fst.csv"),
   "population_1,population_2,FST_Hudson,markers_used",Fst.estimate(d));
  GenotypeData.writeCsv(out.resolve("fis.csv"),
   "population_id,Ho_sum,He_sum,FIS,loci",Fis.estimate(d));
  GenotypeData.writeCsv(out.resolve("ld.csv"),
   "marker_1,marker_2,n_complete,r,r_squared",LinkageDisequilibrium.estimate(d,20));
  GenotypeData.writeCsv(out.resolve("nei_distance.csv"),
   "population_1,population_2,nei_identity,nei_distance,loci",GeneticDistance.estimate(d));
  GenotypeData.writeCsv(out.resolve("amova.csv"),
   "marker_id,n,populations,SS_among,SS_within,var_among,var_within,Phi_among",Amova.estimate(d));
  System.out.println("Completed population-genetics estimators. Haplotype analysis requires phased allele rows.");
 }
}
