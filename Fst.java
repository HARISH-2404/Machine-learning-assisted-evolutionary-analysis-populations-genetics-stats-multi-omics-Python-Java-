package org.ipgml.populationgenetics;
import java.util.*;
/** Pairwise Hudson-style FST ratio of sums; negative estimates are retained. */
public class Fst {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();List<String> pops=d.uniquePopulations();
  for(int a=0;a<pops.size();a++)for(int b=a+1;b<pops.size();b++){
   double num=0,den=0;int used=0;String p1=pops.get(a),p2=pops.get(b);
   for(int m=0;m<d.markers.size();m++){
    double s1=0,s2=0;int n1=0,n2=0;
    for(int i=0;i<d.sampleIds.size();i++){
     double v=d.dosage[i][m];
     if(!Double.isFinite(v))continue;
     if(d.populations.get(i).equals(p1)){s1+=v;n1++;}
     if(d.populations.get(i).equals(p2)){s2+=v;n2++;}
    }
    if(n1<2||n2<2)continue;
    double x=s1/(2*n1),y=s2/(2*n2),dd=x*(1-y)+y*(1-x);
    if(dd<=0)continue;
    num+=(x-y)*(x-y)-x*(1-x)/(n1-1)-y*(1-y)/(n2-1);den+=dd;used++;
   }
   out.add(p1+","+p2+","+(den>0?num/den:Double.NaN)+","+used);
  }return out;
 }
}
