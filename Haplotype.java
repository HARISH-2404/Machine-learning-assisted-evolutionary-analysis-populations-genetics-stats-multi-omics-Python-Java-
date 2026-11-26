package org.ipgml.populationgenetics;
import java.util.*;
/** Haplotype frequencies from phased 0/1 chromosome rows, not unphased dosages. */
public class Haplotype {
 public static List<String> estimate(int[][] phased,String[] markerNames){
  if(phased.length==0||markerNames.length==0)throw new IllegalArgumentException("Empty phased matrix.");
  Map<String,Integer> counts=new TreeMap<>();
  for(int[] row:phased){
   if(row.length!=markerNames.length)throw new IllegalArgumentException("Marker count mismatch.");
   StringBuilder s=new StringBuilder();
   for(int v:row){if(v!=0&&v!=1)throw new IllegalArgumentException("Phased alleles must be 0/1.");s.append(v);}
   counts.put(s.toString(),counts.getOrDefault(s.toString(),0)+1);
  }
  List<String> out=new ArrayList<>();
  for(String h:counts.keySet())out.add(h+","+counts.get(h)+","+(double)counts.get(h)/phased.length);
  return out;
 }
}
