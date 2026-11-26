package org.ipgml.evolutionarygenetics;
import java.util.Random;
/** Wright-Fisher drift. Binomial sampling uses Bernoulli gene-copy draws. */
public final class GeneticDrift {
 private GeneticDrift() {}
 private static int binomial(int n,double p,Random r){int k=0;for(int i=0;i<n;i++)if(r.nextDouble()<p)k++;return k;}
 public static double[][] simulate(double p,int generations,int ne,int replicates,long seed){
  Selection.checkP(p);if(generations<0||ne<1||replicates<1)throw new IllegalArgumentException("Invalid generations/Ne/replicates");
  Random r=new Random(seed);double[][] out=new double[replicates][generations+1];
  for(int j=0;j<replicates;j++){out[j][0]=p;for(int t=1;t<=generations;t++)out[j][t]=binomial(2*ne,out[j][t-1],r)/(2.0*ne);}return out;
 }
}
