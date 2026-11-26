package org.ipgml.genomics;
import java.util.*;
/** Single-marker additive regression GWAS with optional numeric covariates. */
public class Gwas {
  public static class Result { public int marker,n; public double beta,se,t,p,maf,r2; }
  private static double[] solve(double[][] a,double[] b){
    int n=b.length; double[][] m=new double[n][n+1];
    for(int i=0;i<n;i++){System.arraycopy(a[i],0,m[i],0,n);m[i][n]=b[i];}
    for(int k=0;k<n;k++){int pivot=k;for(int i=k+1;i<n;i++)if(Math.abs(m[i][k])>Math.abs(m[pivot][k]))pivot=i;
      if(Math.abs(m[pivot][k])<1e-12)return null;double[] tmp=m[k];m[k]=m[pivot];m[pivot]=tmp;
      double d=m[k][k];for(int j=k;j<=n;j++)m[k][j]/=d;
      for(int i=0;i<n;i++)if(i!=k){double f=m[i][k];for(int j=k;j<=n;j++)m[i][j]-=f*m[k][j];}
    }double[] x=new double[n];for(int i=0;i<n;i++)x[i]=m[i][n];return x;
  }
  private static double[][] inverse(double[][] a){
    int n=a.length;double[][] inv=new double[n][n];
    for(int j=0;j<n;j++){double[] e=new double[n];e[j]=1;double[] x=solve(a,e);if(x==null)return null;for(int i=0;i<n;i++)inv[i][j]=x[i];}
    return inv;
  }
  private static double tTwoSidedApprox(double t,int df){
    // Normal-tail approximation for CLI portability; use a statistical library for exact small-df inference.
    return Math.min(1.0,2.0*(1.0-normalCdf(Math.abs(t))));
  }
  private static double normalCdf(double x){
    double z=Math.abs(x)/Math.sqrt(2),t=1/(1+0.3275911*z);
    double erf=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-z*z);
    return 0.5*(1+Math.copySign(erf,x));
  }
  public static List<Result> scan(double[][] g,double[] y,double[][] cov){
    int n=y.length;if(g.length!=n)throw new IllegalArgumentException("sample dimension mismatch");
    int c=cov==null?0:cov[0].length;List<Result> out=new ArrayList<>();
    for(int j=0;j<g[0].length;j++){
      List<Integer> keep=new ArrayList<>();double alt=0;int called=0,het=0;
      for(int i=0;i<n;i++)if(Double.isFinite(g[i][j])&&Double.isFinite(y[i])){
        boolean ok=true;if(cov!=null)for(double v:cov[i])if(!Double.isFinite(v))ok=false;
        if(ok){keep.add(i);alt+=g[i][j];called++;if(g[i][j]==1)het++;}
      }
      Result r=new Result();r.marker=j;r.n=called;r.maf=called==0?Double.NaN:Math.min(alt/(2*called),1-alt/(2*called));
      int k=c+2;if(called<=k+2){r.beta=r.se=r.t=r.p=r.r2=Double.NaN;out.add(r);continue;}
      double[][] X=new double[called][k];double[] yy=new double[called];int z=0;
      for(int i:keep){X[z][0]=1;if(cov!=null)for(int q=0;q<c;q++)X[z][q+1]=cov[i][q];X[z][k-1]=g[i][j];yy[z]=y[i];z++;}
      double[][] xtx=new double[k][k];double[] xty=new double[k];
      for(int a=0;a<k;a++){for(int b=0;b<k;b++)for(int i=0;i<called;i++)xtx[a][b]+=X[i][a]*X[i][b];for(int i=0;i<called;i++)xty[a]+=X[i][a]*yy[i];}
      double[] beta=solve(xtx,xty);double[][] inv=inverse(xtx);
      if(beta==null||inv==null){r.beta=r.se=r.t=r.p=r.r2=Double.NaN;out.add(r);continue;}
      double rss=0,tss=0,mean=0;for(double v:yy)mean+=v/called;
      for(int i=0;i<called;i++){double fit=0;for(int a=0;a<k;a++)fit+=X[i][a]*beta[a];rss+=(yy[i]-fit)*(yy[i]-fit);tss+=(yy[i]-mean)*(yy[i]-mean);}
      double s2=rss/(called-k);r.beta=beta[k-1];r.se=Math.sqrt(Math.max(0,s2*inv[k-1][k-1]));
      r.t=r.se>0?r.beta/r.se:Double.NaN;r.p=Double.isFinite(r.t)?tTwoSidedApprox(r.t,called-k):Double.NaN;r.r2=tss>0?1-rss/tss:Double.NaN;out.add(r);
    }return out;
  }
}
