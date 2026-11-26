package org.ipgml.populationgenetics;
import java.util.*;
/** Lightweight PCA via power iteration on standardized genotype matrix.
 * Returns sample scores for requested components. Intended for exploration;
 * use a validated numerical library for large-scale production PCA.
 */
public class PopulationStructure {
 public static double[][] pca(GenotypeData d,int components,int iterations){
  if(components<1||iterations<10)throw new IllegalArgumentException("Invalid PCA settings.");
  int n=d.dosage.length,m=d.markers.size();double[][] z=new double[n][m];int kept=0;
  for(int j=0;j<m;j++){
   double sum=0;int count=0;
   for(int i=0;i<n;i++)if(Double.isFinite(d.dosage[i][j])){sum+=d.dosage[i][j];count++;}
   if(count==0)continue;double mean=sum/count,p=mean/2,sd=Math.sqrt(2*p*(1-p));
   if(sd<1e-8)continue;
   for(int i=0;i<n;i++){double v=Double.isFinite(d.dosage[i][j])?d.dosage[i][j]:mean;z[i][kept]=(v-mean)/sd;}
   kept++;
  }
  if(kept==0)throw new IllegalArgumentException("No polymorphic markers for PCA.");
  double[][] scores=new double[n][Math.min(components,n)];
  // Deflate sample covariance matrix after each eigenvector.
  double[][] c=new double[n][n];
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0;for(int k=0;k<kept;k++)s+=z[i][k]*z[j][k];c[i][j]=s;}
  int kmax=scores[0].length;
  for(int pc=0;pc<kmax;pc++){
   double[] v=new double[n];Arrays.fill(v,1.0/Math.sqrt(n));
   for(int it=0;it<iterations;it++){
    double[] nv=new double[n];double norm=0;
    for(int i=0;i<n;i++){for(int j=0;j<n;j++)nv[i]+=c[i][j]*v[j];norm+=nv[i]*nv[i];}
    norm=Math.sqrt(norm);if(norm<1e-12)break;
    for(int i=0;i<n;i++)v[i]=nv[i]/norm;
   }
   double eig=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)eig+=v[i]*c[i][j]*v[j];
   for(int i=0;i<n;i++)scores[i][pc]=v[i]*Math.sqrt(Math.max(eig,0));
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)c[i][j]-=eig*v[i]*v[j];
  }
  return scores;
 }
}
