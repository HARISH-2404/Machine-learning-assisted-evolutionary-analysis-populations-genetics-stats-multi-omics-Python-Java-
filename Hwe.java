package org.ipgml.populationgenetics;
import java.util.*;
/** Pearson Hardy-Weinberg chi-square test, 1 df; asymptotic p via erfc approximation. */
public class Hwe {
 static double erfc(double x){
  double z=Math.abs(x),t=1/(1+0.5*z);
  double ans=t*Math.exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
   t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return x>=0?ans:2-ans;
 }
 public static List<String> test(GenotypeData d){
  List<String> out=new ArrayList<>();
  for(String pop:d.uniquePopulations())for(int m=0;m<d.markers.size();m++){
   double sum=0;int n=0;int[] obs=new int[3];
   for(int i=0;i<d.sampleIds.size();i++)if(d.populations.get(i).equals(pop)){
    double v=d.dosage[i][m];if(Double.isFinite(v)){sum+=v;n++;obs[(int)v]++;}
   }
   if(n==0)continue;
   double p=sum/(2*n),q=1-p;
   double[] e={n*q*q,n*2*p*q,n*p*p}; double chi=0;
   boolean valid=true;
   for(int k=0;k<3;k++){if(e[k]<5)valid=false;if(e[k]>0)chi+=(obs[k]-e[k])*(obs[k]-e[k])/e[k];else chi=Double.NaN;}
   double pv=Double.isFinite(chi)?erfc(Math.sqrt(chi/2)):Double.NaN;
   out.add(pop+","+d.markers.get(m)+","+n+","+p+","+obs[0]+","+obs[1]+","+obs[2]+","+
     e[0]+","+e[1]+","+e[2]+","+chi+","+pv+","+valid);
  }return out;
 }
}
