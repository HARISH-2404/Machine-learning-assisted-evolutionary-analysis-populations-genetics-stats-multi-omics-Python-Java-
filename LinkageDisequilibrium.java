package org.ipgml.populationgenetics;
import java.util.*;
/** Pairwise dosage correlation r and r-squared. Not phased D-prime. */
public class LinkageDisequilibrium {
 public static List<String> estimate(GenotypeData d,int minComplete){
  if(minComplete<3)throw new IllegalArgumentException("minComplete must be >=3.");
  List<String> out=new ArrayList<>();
  for(int a=0;a<d.markers.size();a++)for(int b=a+1;b<d.markers.size();b++){
   double sx=0,sy=0,sxx=0,syy=0,sxy=0;int n=0;
   for(int i=0;i<d.sampleIds.size();i++){
    double x=d.dosage[i][a],y=d.dosage[i][b];
    if(Double.isFinite(x)&&Double.isFinite(y)){n++;sx+=x;sy+=y;sxx+=x*x;syy+=y*y;sxy+=x*y;}
   }
   if(n<minComplete)continue;
   double cov=sxy-sx*sy/n,vx=sxx-sx*sx/n,vy=syy-sy*sy/n;
   if(vx<=0||vy<=0)continue;
   double r=cov/Math.sqrt(vx*vy);
   out.add(d.markers.get(a)+","+d.markers.get(b)+","+n+","+r+","+(r*r));
  }return out;
 }
}
