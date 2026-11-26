package org.ipgml.evolutionarygenetics;
/** Two-way mutation: u A->a, v a->A. */
public final class Mutation {
 private Mutation() {}
 public static double step(double p,double u,double v){Selection.checkP(p);Selection.checkP(u);Selection.checkP(v);return p*(1-u)+(1-p)*v;}
 public static double[][] simulate(double p,int generations,double u,double v){if(generations<0)throw new IllegalArgumentException("generations >=0");double[][] x=new double[generations+1][2];for(int t=0;t<=generations;t++){x[t]=new double[]{t,p};p=step(p,u,v);}return x;}
}
