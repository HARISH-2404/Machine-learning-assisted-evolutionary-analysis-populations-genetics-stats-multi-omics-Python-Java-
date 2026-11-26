package org.ipgml.quantitativegenetics;
import java.util.*;
/** Gaussian-prior Bayesian ridge posterior mean for fixed residual variance and prior variance. */
public final class BayesianRidge {
 public static double[] posteriorMean(double[][] X,double[] y,double residualVariance,double markerVariance){
  if(residualVariance<=0||markerVariance<=0||X.length!=y.length)throw new IllegalArgumentException("Invalid inputs");
  int n=X.length,p=X[0].length;double[][] A=new double[p][p+1];double[] mean=new double[p];for(double[] row:X)for(int j=0;j<p;j++)mean[j]+=row[j]/n;
  double ym=Arrays.stream(y).average().orElseThrow();
  for(int i=0;i<p;i++){for(int j=0;j<p;j++){for(int r=0;r<n;r++)A[i][j]+=(X[r][i]-mean[i])*(X[r][j]-mean[j])/residualVariance;if(i==j)A[i][j]+=1/markerVariance;}
   for(int r=0;r<n;r++)A[i][p]+=(X[r][i]-mean[i])*(y[r]-ym)/residualVariance;}
  for(int k=0;k<p;k++){double d=A[k][k];if(Math.abs(d)<1e-12)throw new IllegalArgumentException("Singular posterior precision");for(int j=k;j<=p;j++)A[k][j]/=d;for(int i=0;i<p;i++)if(i!=k){d=A[i][k];for(int j=k;j<=p;j++)A[i][j]-=d*A[k][j];}}
  double[] b=new double[p];for(int j=0;j<p;j++)b[j]=A[j][p];return b;
 }
}
