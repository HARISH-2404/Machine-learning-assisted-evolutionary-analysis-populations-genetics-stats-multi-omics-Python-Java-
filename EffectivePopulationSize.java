package org.ipgml.evolutionarygenetics;
/** Temporal Ne approximation and neutral heterozygosity decay. */
public final class EffectivePopulationSize {
 private EffectivePopulationSize() {}
 public static double estimate(double p0,double p1,int generations){if(!(p0>0&&p0<1&&p1>=0&&p1<=1)||generations<1)throw new IllegalArgumentException("Invalid frequencies/generations");double f=(p1-p0)*(p1-p0)/(p0*(1-p0));return f==0?Double.POSITIVE_INFINITY:generations/(2*f);}
 public static double heterozygosityAfter(double h0,int ne,int generations){if(h0<0||h0>1||ne<1||generations<0)throw new IllegalArgumentException("Invalid H0/Ne/generations");return h0*Math.pow(1-1.0/(2*ne),generations);}
}
