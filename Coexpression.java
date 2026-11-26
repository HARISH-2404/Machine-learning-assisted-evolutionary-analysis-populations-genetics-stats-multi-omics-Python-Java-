package org.ipgml.transcriptomics;
import java.util.*;
/** Pairwise Pearson co-expression network edge extraction. */
public class Coexpression {
  public static class Edge {public int gene1,gene2;public double correlation;}
  public static List<Edge> edges(double[][] expression,double threshold,boolean absolute){
    int n=expression.length,p=expression[0].length;List<Edge> out=new ArrayList<>();
    for(int a=0;a<p;a++)for(int b=a+1;b<p;b++){
      double sx=0,sy=0,sxx=0,syy=0,sxy=0;int c=0;
      for(int i=0;i<n;i++){double x=expression[i][a],y=expression[i][b];if(Double.isFinite(x)&&Double.isFinite(y)){sx+=x;sy+=y;sxx+=x*x;syy+=y*y;sxy+=x*y;c++;}}
      double dx=c*sxx-sx*sx,dy=c*syy-sy*sy;if(c<4||dx<=0||dy<=0)continue;
      double r=(c*sxy-sx*sy)/Math.sqrt(dx*dy);
      if(absolute?Math.abs(r)>=threshold:r>=threshold){Edge e=new Edge();e.gene1=a;e.gene2=b;e.correlation=r;out.add(e);}
    }return out;
  }
}
