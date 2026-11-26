package org.ipgml.integration;
import java.util.*;
public final class IntegrationUtils {
 private IntegrationUtils(){}
 public static List<String> intersectSamples(List<OmicsMatrix> blocks){
  if(blocks.isEmpty()) return List.of(); Set<String> s=new LinkedHashSet<>(blocks.get(0).samples);
  for(OmicsMatrix b:blocks) s.retainAll(b.samples); return new ArrayList<>(s);
 }
 public static double[][] zScore(double[][] train,double[][] test){
  if(train.length==0) throw new IllegalArgumentException("Empty training matrix"); int p=train[0].length;
  double[] mean=new double[p],sd=new double[p];
  for(double[] r:train) for(int j=0;j<p;j++) mean[j]+=r[j]/train.length;
  for(double[] r:train) for(int j=0;j<p;j++) sd[j]+=Math.pow(r[j]-mean[j],2)/train.length;
  for(int j=0;j<p;j++) sd[j]=Math.sqrt(sd[j])==0?1:Math.sqrt(sd[j]);
  double[][] out=new double[test.length][p];
  for(int i=0;i<test.length;i++){if(test[i].length!=p) throw new IllegalArgumentException("Feature mismatch"); for(int j=0;j<p;j++) out[i][j]=(test[i][j]-mean[j])/sd[j];}
  return out;
 }
 public static double[] weightedLateIntegration(double[][] scores,double[] weights){
  if(scores.length!=weights.length || scores.length==0) throw new IllegalArgumentException("One weight per omics score block required");
  double sum=0; for(double w:weights){if(w<0) throw new IllegalArgumentException("Negative weight"); sum+=w;} if(sum==0) throw new IllegalArgumentException("Weights sum to zero");
  int n=scores[0].length; double[] out=new double[n];
  for(int k=0;k<scores.length;k++){if(scores[k].length!=n) throw new IllegalArgumentException("Sample order/length mismatch"); for(int i=0;i<n;i++) out[i]+=scores[k][i]*weights[k]/sum;}
  return out;
 }
}
