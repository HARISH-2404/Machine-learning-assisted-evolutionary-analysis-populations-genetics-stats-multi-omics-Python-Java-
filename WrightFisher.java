package org.ipgml.forecasting;
import java.util.*;
/** Diploid single-locus Wright-Fisher simulator. */
public final class WrightFisher {
 private WrightFisher(){}
 public static double[] fitness(double s,double h){return new double[]{1+s,1+h*s,1};}
 public static double select(double p,double[] w){
  double q=1-p,mean=p*p*w[0]+2*p*q*w[1]+q*q*w[2];
  if(p<0||p>1||mean<=0) throw new IllegalArgumentException("Invalid p or fitness");
  return (p*p*w[0]+p*q*w[1])/mean;
 }
 public static double[][] simulate(double p0,int generations,int populationSize,double s,double h,int replicates,long seed){
  if(generations<0||populationSize<1||replicates<1) throw new IllegalArgumentException("Invalid simulation dimensions");
  Random rng=new Random(seed); double[][] out=new double[replicates][generations+1];
  for(int r=0;r<replicates;r++){double p=p0;out[r][0]=p;
   for(int g=1;g<=generations;g++){p=select(p,fitness(s,h));int copies=0;
    for(int k=0;k<2*populationSize;k++) if(rng.nextDouble()<p) copies++;
    p=(double)copies/(2*populationSize);out[r][g]=p;
   }
  } return out;
 }
}
