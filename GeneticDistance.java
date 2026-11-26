package org.ipgml.populationgenetics;
import java.util.*;
/** Nei standard distance D=-ln(I), using both alleles at each biallelic locus. */
public class GeneticDistance {
 public static List<String> estimate(GenotypeData d){
  List<String> out=new ArrayList<>();List<String> pops=d.uniquePopulations();
  for(int a=0;a<pops.size();a++)for(int b=a+1;b<pops.size();b++){
   String p1=pops.get(a),p2=pops.get(b);double dot=0,n1=0,n2=0;int loci=0;
   for(int m=0;m<d.markers.size();m++){
    double s1=0,s2=0;int c1=0,c2=0;
    for(int i=0;i<d.sampleIds.size();i++){double v=d.dosage[i][m];if(!Double.isFinite(v))continue;
     if(d.populations.get(i).equals(p1)){s1+=v;c1++;}
     if(d.populations.get(i).equals(p2)){s2+=v;c2++;}}
    if(c1==0||c2==0)continue;
    double x=s1/(2*c1),y=s2/(2*c2);
    double[] aa={x,1-x},bb={y,1-y};
    for(int k=0;k<2;k++){dot+=Math.sqrt(aa[k]*bb[k]);n1+=aa[k]*aa[k];n2+=bb[k]*bb[k];}
    loci++;
   }
   double identity=dot/Math.sqrt(n1*n2),dist=-Math.log(Math.max(identity,1e-15));
   out.add(p1+","+p2+","+identity+","+dist+","+loci);
  }return out;
 }
}
