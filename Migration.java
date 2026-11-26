package org.ipgml.evolutionarygenetics;
/** Deterministic migration; matrix rows=source and columns=destination. */
public final class Migration {
 private Migration() {}
 public static double[] step(double[] p,double[][] m){int n=p.length;if(m.length!=n)throw new IllegalArgumentException("matrix dimensions");for(double x:p)Selection.checkP(x);double[] out=new double[n];for(int i=0;i<n;i++){if(m[i].length!=n)throw new IllegalArgumentException("matrix must be square");double sum=0;for(int j=0;j<n;j++){if(m[i][j]<0||m[i][j]>1)throw new IllegalArgumentException("migration entries in [0,1]");sum+=m[i][j];out[j]+=p[i]*m[i][j];}if(Math.abs(sum-1)>1e-9)throw new IllegalArgumentException("each row must sum to one");}return out;}
 public static double[][] simulate(double[] p,double[][] m,int generations){if(generations<0)throw new IllegalArgumentException("generations >=0");double[][] out=new double[generations+1][p.length];for(int t=0;t<=generations;t++){out[t]=p.clone();if(t<generations)p=step(p,m);}return out;}
}
