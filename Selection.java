package org.ipgml.evolutionarygenetics;
/** Deterministic viability selection. */
public final class Selection {
 private Selection() {}
 public static double[] step(double p,double[] w) {
  checkP(p); if(w==null||w.length!=3) throw new IllegalArgumentException("fitness=[wAA,wAa,waa]");
  double q=1-p, aa=p*p, ab=2*p*q, bb=q*q;
  double mean=aa*w[0]+ab*w[1]+bb*w[2]; if(mean<=0) throw new IllegalArgumentException("Mean fitness is zero");
  return new double[]{(aa*w[0]+0.5*ab*w[1])/mean,mean};
 }
 public static double[][] simulate(double p,int generations,double[] fitness) {
  if(generations<0) throw new IllegalArgumentException("generations must be >=0");
  double[][] out=new double[generations+1][3];
  for(int t=0;t<=generations;t++){ double[] z=step(p,fitness); out[t]=new double[]{t,p,z[1]}; p=z[0]; }
  return out;
 }
 static void checkP(double p){if(!Double.isFinite(p)||p<0||p>1)throw new IllegalArgumentException("p must be in [0,1]");}
}
