package org.ipgml.evolutionarygenetics;
import java.util.Random;
/** Selection -> mutation -> optional Wright-Fisher sampling each generation. */
public final class AlleleTrajectories {
 private AlleleTrajectories() {}
 public static double[][] simulate(double p,int generations,double[] fitness,double u,double v,int ne,long seed){
  Selection.checkP(p);if(generations<0||ne<0)throw new IllegalArgumentException("generations >=0; Ne=0 means deterministic");
  Random r=new Random(seed);double[][] out=new double[generations+1][2];
  for(int t=0;t<=generations;t++){out[t]=new double[]{t,p};if(t<generations){p=Selection.step(p,fitness)[0];p=Mutation.step(p,u,v);if(ne>0){int k=0;for(int i=0;i<2*ne;i++)if(r.nextDouble()<p)k++;p=k/(2.0*ne);}}}return out;
 }
}
