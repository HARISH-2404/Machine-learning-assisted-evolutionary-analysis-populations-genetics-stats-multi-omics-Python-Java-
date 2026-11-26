package org.ipgml.integration;
import java.util.*;
public class IntegrationCli {
 public static void main(String[] args){
  double[][] train={{1,2},{3,4}}, test={{5,6}};
  double[][] scaled=IntegrationUtils.zScore(train,test);
  System.out.printf("Training-statistics z-score: %.3f, %.3f%n",scaled[0][0],scaled[0][1]);
  List<LocusObservation> series=List.of(
   new LocusObservation("Chr1",100,"GeneA","A","Pop1","Site1","T1",1,2,10),
   new LocusObservation("Chr1",100,"GeneA","A","Pop1","Site1","T2",3,8,10));
  System.out.println("chromosome\\tposition\\tgene\\tallele\\tpopulation\\tlocation\\ttime\\tgeneration\\tfreq\\tdelta\\tdelta_per_generation");
  EvolutionaryTrajectory.summarize(series).forEach(System.out::println);
 }
}
