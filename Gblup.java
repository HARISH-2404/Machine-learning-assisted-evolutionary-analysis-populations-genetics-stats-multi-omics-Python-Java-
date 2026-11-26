package org.ipgml.quantitativegenetics;
import java.util.*;
/** Kernel GBLUP with supplied GRM and fixed lambda; no REML estimation. */
public final class Gblup {
 public static Map<String,Object> predict(double[][] G,double[] y,double lambda){
  int n=y.length;if(G.length!=n||lambda<=0)throw new IllegalArgumentException("Invalid G/y/lambda");
  double[][] A=new double[n][n+1];for(int i=0;i<n;i++){if(G[i].length!=n)throw new IllegalArgumentException("G must be square");for(int j=0;j<n;j++)A[i][j]=G[i][j]+(i==j?lambda:0);A[i][n]=y[i];}
  // Gauss-Jordan solve A z = y
  for(int k=0;k<n;k++){int pivot=k;for(int i=k+1;i<n;i++)if(Math.abs(A[i][k])>Math.abs(A[pivot][k]))pivot=i;double[]tmp=A[k];A[k]=A[pivot];A[pivot]=tmp;if(Math.abs(A[k][k])<1e-12)throw new IllegalArgumentException("Singular G+lambda I");double d=A[k][k];for(int j=k;j<=n;j++)A[k][j]/=d;for(int i=0;i<n;i++)if(i!=k){d=A[i][k];for(int j=k;j<=n;j++)A[i][j]-=d*A[k][j];}}
  double[] alpha=new double[n],u=new double[n],pred=new double[n];double mean=Arrays.stream(y).average().orElseThrow();for(int i=0;i<n;i++)alpha[i]=A[i][n];for(int i=0;i<n;i++){for(int j=0;j<n;j++)u[i]+=G[i][j]*alpha[j];pred[i]=mean+u[i];}
  Map<String,Object>out=new LinkedHashMap<>();out.put("genomic_values",u);out.put("predictions",pred);out.put("lambda",lambda);return out;
 }
}
