package org.ipgml.evolutionarygenetics;
/** Basic SFS summaries. Input rows are [reference count, alternate count]. */
public final class SelectionSweep {
 private SelectionSweep() {}
 public static double nucleotideDiversity(double[][] counts){if(counts.length==0)return Double.NaN;double sum=0;for(double[] c:counts){if(c.length!=2||c[0]<0||c[1]<0)throw new IllegalArgumentException("Each site needs two nonnegative counts");double n=c[0]+c[1];if(n>1)sum+=2*c[0]*c[1]/(n*(n-1));}return sum/counts.length;}
 public static double tajimasD(double[][] c){if(c.length<2) return Double.NaN;int n=(int)(c[0][0]+c[0][1]);if(n<2)throw new IllegalArgumentException("sample size >=2");int S=0;double pi=0;for(double[] x:c){if(x.length!=2||(int)(x[0]+x[1])!=n)throw new IllegalArgumentException("Equal sample size required");if(x[0]>0&&x[1]>0)S++;pi+=2*x[0]*x[1]/(n*(n-1));}if(S<2)return Double.NaN;double a1=0,a2=0;for(int i=1;i<n;i++){a1+=1.0/i;a2+=1.0/(i*i);}double b1=(n+1.0)/(3*(n-1)),b2=2.0*(n*n+n+3)/(9*n*(n-1));double c1=b1-1/a1,c2=b2-(n+2.0)/(a1*n)+a2/(a1*a1);double e1=c1/a1,e2=c2/(a1*a1+a2);double den=Math.sqrt(e1*S+e2*S*(S-1));return den==0?Double.NaN:(pi-S/a1)/den;}
}
