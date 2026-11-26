package org.ipgml.integration;
/** Transparent baseline population-genetic estimators; diploid biallelic assumptions. */
public final class EvolutionaryEstimators {
 private EvolutionaryEstimators(){}
 public static double alleleFrequency(long alleleCount,long totalAlleles){
  if(totalAlleles<=0 || alleleCount<0 || alleleCount>totalAlleles) throw new IllegalArgumentException("Invalid allele counts");
  return (double)alleleCount/totalAlleles;
 }
 public static double expectedHeterozygosity(double p){return 2*p*(1-p);}
 public static double deterministicSelectionStep(double p,double s,double h){
  if(p<0||p>1) throw new IllegalArgumentException("p outside [0,1]");
  double q=1-p,wAA=1+s,wAa=1+h*s,waa=1;
  if(Math.min(wAA,Math.min(wAa,waa))<0) throw new IllegalArgumentException("Negative fitness");
  double mean=p*p*wAA+2*p*q*wAa+q*q*waa;
  return (p*p*wAA+p*q*wAa)/mean;
 }
 public static double temporalNeRough(double p0,double p1,double generations){
  double den=2*Math.pow(p1-p0,2);
  return den==0?Double.NaN:p0*(1-p0)*generations/den;
 }
 public static double fstApprox(double p1,double p2){
  double mean=(p1+p2)/2,ht=2*mean*(1-mean);
  double hs=(2*p1*(1-p1)+2*p2*(1-p2))/2;
  return ht==0?Double.NaN:(ht-hs)/ht;
 }
 public static double selectionLogitSlope(double[] generations,double[] frequencies){
  if(generations.length!=frequencies.length||generations.length<2) throw new IllegalArgumentException("Need >=2 matched time points");
  double mx=0,my=0,eps=1e-6;
  for(int i=0;i<generations.length;i++){double p=Math.max(eps,Math.min(1-eps,frequencies[i])); mx+=generations[i]/generations.length; my+=Math.log(p/(1-p))/generations.length;}
  double cov=0,var=0;
  for(int i=0;i<generations.length;i++){double p=Math.max(eps,Math.min(1-eps,frequencies[i])); double y=Math.log(p/(1-p)); cov+=(generations[i]-mx)*(y-my); var+=(generations[i]-mx)*(generations[i]-mx);}
  return var==0?Double.NaN:cov/var;
 }
}
