package org.ipgml.genomics;
import java.util.*;
/** Marker-grid additive QTL scan; not composite interval mapping. */
public class Qtl {
  public static class Result {public int marker,n;public double beta, lod, fStatistic;}
  public static List<Result> scan(double[][] g,double[] y){
    List<Result> out=new ArrayList<>();int n=y.length;
    for(int j=0;j<g[0].length;j++){
      double sy=0,sx=0,sxx=0,sxy=0;int count=0;
      for(int i=0;i<n;i++)if(Double.isFinite(g[i][j])&&Double.isFinite(y[i])){double x=g[i][j];sx+=x;sy+=y[i];sxx+=x*x;sxy+=x*y[i];count++;}
      Result r=new Result();r.marker=j;r.n=count;
      double den=count*sxx-sx*sx;
      if(count<6||den<=0){r.beta=r.lod=r.fStatistic=Double.NaN;out.add(r);continue;}
      double beta=(count*sxy-sx*sy)/den,intercept=(sy-beta*sx)/count;
      double mean=sy/count,rss1=0,rss0=0;
      for(int i=0;i<n;i++)if(Double.isFinite(g[i][j])&&Double.isFinite(y[i])){
        double e=y[i]-intercept-beta*g[i][j];rss1+=e*e;rss0+=(y[i]-mean)*(y[i]-mean);
      }
      r.beta=beta;r.fStatistic=((rss0-rss1)/(rss1/(count-2)));
      r.lod=(rss1>0&&rss0>0)?count/2.0*Math.log10(rss0/rss1):Double.NaN;out.add(r);
    }return out;
  }
}
