package org.ipgml.machinelearning;
import java.util.*;
public final class Classification {
 public static double sigmoid(double z){return 1/(1+Math.exp(-Math.max(-35,Math.min(35,z))));}
 public static Map<String,Object> logistic(double[][] X,int[] y,double rate,int epochs,double l2){
  int n=X.length;if(n==0||y.length!=n)throw new IllegalArgumentException("Invalid inputs");int p=X[0].length;double[] mu=new double[p],sd=new double[p],w=new double[p+1];
  for(double[] row:X)for(int j=0;j<p;j++)mu[j]+=row[j]/n;for(double[] row:X)for(int j=0;j<p;j++)sd[j]+=(row[j]-mu[j])*(row[j]-mu[j])/n;for(int j=0;j<p;j++){sd[j]=Math.sqrt(sd[j]);if(sd[j]==0)sd[j]=1;}
  for(int it=0;it<epochs;it++){double[] grad=new double[p+1];for(int i=0;i<n;i++){double z=w[0];for(int j=0;j<p;j++)z+=w[j+1]*(X[i][j]-mu[j])/sd[j];double err=sigmoid(z)-y[i];grad[0]+=err;for(int j=0;j<p;j++)grad[j+1]+=err*(X[i][j]-mu[j])/sd[j];}for(int j=0;j<=p;j++){grad[j]/=n;if(j>0)grad[j]+=l2*w[j];w[j]-=rate*grad[j];}}
  double[] prob=new double[n];int[] pred=new int[n];for(int i=0;i<n;i++){double z=w[0];for(int j=0;j<p;j++)z+=w[j+1]*(X[i][j]-mu[j])/sd[j];prob[i]=sigmoid(z);pred[i]=prob[i]>=.5?1:0;}
  Map<String,Object>out=new LinkedHashMap<>();out.put("coefficients",Arrays.copyOfRange(w,1,w.length));out.put("intercept",w[0]);out.put("probabilities",prob);out.put("predictions",pred);return out;
 }
}
