package org.ipgml.transcriptomics;
import java.io.*;
import java.util.*;

/** Simple two-group Welch test on log2(CPM+prior); exploratory only. */
public class DifferentialExpression {
  private static double normalCdf(double x){
    double z=Math.abs(x)/Math.sqrt(2),t=1/(1+0.3275911*z);
    double erf=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-z*z);
    return 0.5*(1+Math.copySign(erf,x));
  }
  private static double mean(double[] x,int[] idx){double s=0;for(int i:idx)s+=x[i];return s/idx.length;}
  private static double variance(double[] x,int[] idx,double m){double s=0;for(int i:idx)s+=(x[i]-m)*(x[i]-m);return s/(idx.length-1);}
  public static void run(ExpressionMatrix m,String[] groups,String control,String treatment,String output)throws IOException{
    if(groups.length!=m.samples.size())throw new IllegalArgumentException("groups must match sample columns");
    List<Integer> c=new ArrayList<>(),t=new ArrayList<>();
    for(int i=0;i<groups.length;i++){if(groups[i].equals(control))c.add(i);if(groups[i].equals(treatment))t.add(i);}
    if(c.size()<2||t.size()<2)throw new IllegalArgumentException("At least two biological replicates per group");
    int[] ci=c.stream().mapToInt(Integer::intValue).toArray(),ti=t.stream().mapToInt(Integer::intValue).toArray();
    double[][] log=m.cpm(true,0.5);List<String> lines=new ArrayList<>();
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      w.write("gene_id,mean_control,mean_treatment,log2_fold_change,t_statistic,p_value\n");
      for(int g=0;g<m.genes.size();g++){double[] row=log[g];double mc=mean(row,ci),mt=mean(row,ti);
        double vc=variance(row,ci,mc),vt=variance(row,ti,mt);
        double se=Math.sqrt(vc/ci.length+vt/ti.length);double stat=se>0?(mt-mc)/se:Double.NaN;
        double p=Double.isFinite(stat)?2*(1-normalCdf(Math.abs(stat))):Double.NaN;
        w.write(m.genes.get(g)+","+mc+","+mt+","+(mt-mc)+","+stat+","+p+"\n");
      }
    }
  }
}
