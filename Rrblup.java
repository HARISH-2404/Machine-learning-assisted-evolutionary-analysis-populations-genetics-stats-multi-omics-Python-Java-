package org.ipgml.quantitativegenetics;
import java.util.*;
/** RR-BLUP marker regression; standardized X and fixed ridge penalty. */
public final class Rrblup {
 public static Map<String,Object> fit(double[][] X,double[] y,double lambda){
  int n=X.length;if(n!=y.length||n==0||lambda<=0)throw new IllegalArgumentException("Invalid dimensions/lambda");int p=X[0].length;
  double[] mean=new double[p],sd=new double[p];for(double[] row:X){if(row.length!=p)throw new IllegalArgumentException("Ragged X");for(int j=0;j<p;j++)mean[j]+=row[j]/n;}
  for(double[] row:X)for(int j=0;j<p;j++)sd[j]+=(row[j]-mean[j])*(row[j]-mean[j])/n;for(int j=0;j<p;j++){sd[j]=Math.sqrt(sd[j]);if(sd[j]==0)sd[j]=1;}
  double[][] A=new double[p][p+1];double ym=Arrays.stream(y).average().orElseThrow();
  for(int i=0;i<p;i++){for(int j=0;j<p;j++){for(int r=0;r<n;r++)A[i][j]+=(X[r][i]-mean[i])/sd[i]*(X[r][j]-mean[j])/sd[j];if(i==j)A[i][j]+=lambda;}
   for(int r=0;r<n;r++)A[i][p]+=(X[r][i]-mean[i])/sd[i]*(y[r]-ym);}
  for(int k=0;k<p;k++){double d=A[k][k];if(Math.abs(d)<1e-12)throw new IllegalArgumentException("Singular system");for(int j=k;j<=p;j++)A[k][j]/=d;for(int i=0;i<p;i++)if(i!=k){d=A[i][k];for(int j=k;j<=p;j++)A[i][j]-=d*A[k][j];}}
  double[] b=new double[p],pred=new double[n];for(int j=0;j<p;j++)b[j]=A[j][p];for(int i=0;i<n;i++){pred[i]=ym;for(int j=0;j<p;j++)pred[i]+=(X[i][j]-mean[j])/sd[j]*b[j];}
  Map<String,Object>out=new LinkedHashMap<>();out.put("marker_effects",b);out.put("predictions",pred);out.put("intercept",ym);return out;
 }
}
