package org.ipgml.populationgenetics;
import java.util.*;
/** Counts dosage classes 0,1,2 and estimates their observed frequencies. */
public class GenotypeFrequency {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();
  for(String pop:d.uniquePopulations())for(int m=0;m<d.markers.size();m++){
   int[] c=new int[3];int n=0;
   for(int i=0;i<d.sampleIds.size();i++)if(d.populations.get(i).equals(pop)){
    double v=d.dosage[i][m];if(Double.isFinite(v)){c[(int)v]++;n++;}
   }
   out.add(pop+","+d.markers.get(m)+","+n+","+c[0]+","+c[1]+","+c[2]+","+
    GenotypeData.fmt(n>0?(double)c[0]/n:Double.NaN)+","+
    GenotypeData.fmt(n>0?(double)c[1]/n:Double.NaN)+","+
    GenotypeData.fmt(n>0?(double)c[2]/n:Double.NaN));
  }return out;
 }
}
