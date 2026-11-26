package org.ipgml.integration;
import java.util.*;
/** Calculates simple sequential allele-frequency changes for one ordered series. */
public final class EvolutionaryTrajectory {
 private EvolutionaryTrajectory(){}
 public static List<String> summarize(List<LocusObservation> observations){
  List<LocusObservation> x=new ArrayList<>(observations);
  x.sort(Comparator.comparingLong(o->o.generation));
  List<String> rows=new ArrayList<>(); double prev=Double.NaN; long prevGen=-1;
  for(LocusObservation o:x){double f=o.frequency(); String change=Double.isNaN(prev)?"NA":String.format(Locale.ROOT,"%.6f",f-prev);
   String perGen=(Double.isNaN(prev)||o.generation==prevGen)?"NA":String.format(Locale.ROOT,"%.6f",(f-prev)/(o.generation-prevGen));
   rows.add(String.join("\\t",o.chromosome,Long.toString(o.position),o.geneId,o.allele,o.population,o.location,o.timepoint,Long.toString(o.generation),String.format(Locale.ROOT,"%.6f",f),change,perGen));
   prev=f; prevGen=o.generation;
  } return rows;
 }
}
