package org.ipgml.populationgenetics;
import java.util.*;
/** One-level dosage-based ANOVA variance partition (AMOVA approximation). */
public class Amova {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();List<String> pops=d.uniquePopulations();
  for(int m=0;m<d.markers.size();m++){
   double grandSum=0;int n=0,k=0;Map<String,List<Double>> vals=new LinkedHashMap<>();
   for(String p:pops){List<Double> a=new ArrayList<>();
    for(int i=0;i<d.sampleIds.size();i++)if(d.populations.get(i).equals(p)&&Double.isFinite(d.dosage[i][m]))a.add(d.dosage[i][m]);
    if(!a.isEmpty()){vals.put(p,a);for(double v:a){grandSum+=v;n++;}k++;}
   }
   if(k<2||n<=k)continue;double grand=grandSum/n,ssb=0,ssw=0,sumN2=0;
   for(List<Double> a:vals.values()){double mean=a.stream().mapToDouble(v->v).average().orElse(0);
    ssb+=a.size()*Math.pow(mean-grand,2);sumN2+=(double)a.size()*a.size();
    for(double v:a)ssw+=Math.pow(v-mean,2);}
   int dfb=k-1,dfw=n-k;double msb=ssb/dfb,msw=ssw/dfw;
   double n0=(n-sumN2/n)/dfb,va=(msb-msw)/n0,vw=msw;
   double phi=(va+vw)>0?va/(va+vw):Double.NaN;
   out.add(d.markers.get(m)+","+n+","+k+","+ssb+","+ssw+","+va+","+vw+","+phi);
  }return out;
 }
}
