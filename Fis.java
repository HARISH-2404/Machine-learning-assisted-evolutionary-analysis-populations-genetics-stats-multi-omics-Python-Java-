package org.ipgml.populationgenetics;
import java.util.*;
/** FIS = 1 - sum(Ho)/sum(He), calculated within each population. */
public class Fis {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();
  for(String pop:d.uniquePopulations()){
   double hoSum=0,heSum=0;int loci=0;
   for(int m=0;m<d.markers.size();m++){
    double sum=0;int n=0,het=0;
    for(int i=0;i<d.sampleIds.size();i++)if(d.populations.get(i).equals(pop)){
     double v=d.dosage[i][m];if(Double.isFinite(v)){sum+=v;n++;if(v==1)het++;}
    }
    if(n==0)continue;double p=sum/(2*n),he=2*p*(1-p),ho=(double)het/n;
    if(he>0){hoSum+=ho;heSum+=he;loci++;}
   }
   out.add(pop+","+hoSum+","+heSum+","+(heSum>0?1-hoSum/heSum:Double.NaN)+","+loci);
  }return out;
 }
}
