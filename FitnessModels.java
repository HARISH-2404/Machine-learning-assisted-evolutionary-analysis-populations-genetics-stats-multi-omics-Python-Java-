package org.ipgml.evolutionarygenetics;
/** Diploid fitness triplets ordered [wAA,wAa,waa]. */
public final class FitnessModels {
 private FitnessModels() {}
 public static double[] genotypeFitness(double s,double h) {
  if(s < -1 || h<0 || h>1) throw new IllegalArgumentException("Require s>=-1 and h in [0,1]");
  return new double[]{1+s,1+h*s,1};
 }
 public static double[] overdominance(double wAA,double wAa,double waa) {
  validate(wAA,wAa,waa); if(!(wAa>wAA && wAa>waa)) throw new IllegalArgumentException("Require wAa greater than both homozygotes");
  return new double[]{wAA,wAa,waa};
 }
 public static double[] underdominance(double wAA,double wAa,double waa) {
  validate(wAA,wAa,waa); if(!(wAa<wAA && wAa<waa)) throw new IllegalArgumentException("Require wAa lower than both homozygotes");
  return new double[]{wAA,wAa,waa};
 }
 private static void validate(double a,double b,double c) { if(a<0||b<0||c<0||!Double.isFinite(a+b+c)) throw new IllegalArgumentException("Fitnesses must be finite and nonnegative"); }
}
