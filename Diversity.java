package org.ipgml.populationgenetics;
import java.util.*;
/** Gene diversity and effective allele number Ne=1/(1-He) for biallelic loci. */
public class Diversity {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();
  for(String pop:d.uniquePopulations())for(int m=0;m<d.markers.size();m++){
   double sum=0;int n=0,het=0;
   for(int i=0;i<d.sampleIds.size();i++)if(d.populations.get(i).equals(pop)){
    double v=d.dosage[i][m];if(Double.isFinite(v)){sum+=v;n++;if(v==1)het++;}
   }
   double p=n>0?sum/(2*n):Double.NaN,he=2*p*(1-p),ho=n>0?(double)het/n:Double.NaN;
   double ne=Double.isFinite(he)&&he<1?1/(1-he):Double.NaN;
   out.add(pop+","+d.markers.get(m)+","+ho+","+he+","+ne);
  }return out;
 }
}
