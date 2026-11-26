package org.ipgml.populationgenetics;
import java.util.*;
/** Biallelic Botstein PIC = 2pq(1-pq). */
public class Pic {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();
  for(String pop:d.uniquePopulations())for(int m=0;m<d.markers.size();m++){
   double sum=0;int n=0;
   for(int i=0;i<d.sampleIds.size();i++)if(d.populations.get(i).equals(pop)){
    double v=d.dosage[i][m];if(Double.isFinite(v)){sum+=v;n++;}
   }
   double p=n>0?sum/(2*n):Double.NaN,q=1-p,pic=2*p*q*(1-p*q);
   out.add(pop+","+d.markers.get(m)+","+n+","+p+","+q+","+pic);
  }return out;
 }
}
