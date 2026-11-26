package org.ipgml.machinelearning;
import java.util.*;
public final class Regression {
 public static Map<String,Object> ridge(double[][] X,double[] y,double alpha){
  int n=X.length;if(n==0||y.length!=n||alpha<0)throw new IllegalArgumentException("Invalid inputs");int p=X[0].length;
  double[] mu=new double[p],sd=new double[p];for(double[] row:X)for(int j=0;j<p;j++)mu[j]+=row[j]/n;
  for(double[] row:X)for(int j=0;j<p;j++)sd[j]+=(row[j]-mu[j])*(row[j]-mu[j])/n;for(int j=0;j<p;j++){sd[j]=Math.sqrt(sd[j]);if(sd[j]==0)sd[j]=1;}
  double ym=Arrays.stream(y).average().orElseThrow();double[][] A=new double[p][p+1];
  for(int i=0;i<p;i++){for(int j=0;j<p;j++){for(int r=0;r<n;r++)A[i][j]+=(X[r][i]-mu[i])/sd[i]*(X[r][j]-mu[j])/sd[j];if(i==j)A[i][j]+=alpha;}
   for(int r=0;r<n;r++)A[i][p]+=(X[r][i]-mu[i])/sd[i]*(y[r]-ym);}
  for(int k=0;k<p;k++){double d=A[k][k];if(Math.abs(d)<1e-12)throw new IllegalArgumentException("Singular ridge system");for(int j=k;j<=p;j++)A[k][j]/=d;for(int i=0;i<p;i++)if(i!=k){d=A[i][k];for(int j=k;j<=p;j++)A[i][j]-=d*A[k][j];}}
  double[] b=new double[p],pred=new double[n];for(int j=0;j<p;j++)b[j]=A[j][p];for(int i=0;i<n;i++){pred[i]=ym;for(int j=0;j<p;j++)pred[i]+=(X[i][j]-mu[j])/sd[j]*b[j];}
  Map<String,Object> out=new LinkedHashMap<>();out.put("coefficients",b);out.put("intercept",ym);out.put("predictions",pred);return out;
 }
}
