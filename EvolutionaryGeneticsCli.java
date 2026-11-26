package org.ipgml.evolutionarygenetics;
import java.io.*;
/** Minimal reproducible CLI: java ...EvolutionaryGeneticsCli [output.csv] */
public class EvolutionaryGeneticsCli {
 public static void main(String[] args) throws Exception {
  String file=args.length>0?args[0]:"evolutionary_trajectory.csv";
  double[] w=FitnessModels.genotypeFitness(0.08,0.5);
  double[][] trajectory=AlleleTrajectories.simulate(0.10,100,w,1e-5,1e-6,500,20261127L);
  try(PrintWriter out=new PrintWriter(new FileWriter(file))){out.println("generation,allele_frequency_A");for(double[] row:trajectory)out.printf(java.util.Locale.US,"%.0f,%.8f%n",row[0],row[1]);}
  System.out.println("Wrote "+file+" ("+trajectory.length+" generations)");
 }
}
