package org.ipgml.transcriptomics;
import java.util.*;
/** RNA-seq library QC and normalization helpers. */
public class RnaSeqAnalysis {
  public static double[][] cpm(ExpressionMatrix m){
    double[] lib=m.librarySizes();double[][] out=new double[m.counts.length][m.samples.size()];
    for(int i=0;i<m.counts.length;i++)for(int j=0;j<m.samples.size();j++){
      if(lib[j]<=0)throw new IllegalArgumentException("zero library size");
      out[i][j]=m.counts[i][j]*1e6/lib[j];
    }return out;
  }
  public static int[] detectedGenes(ExpressionMatrix m){
    int[] out=new int[m.samples.size()];
    for(int j=0;j<out.length;j++)for(double[] row:m.counts)if(row[j]>0)out[j]++;
    return out;
  }
}
