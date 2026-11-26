package org.ipgml.transcriptomics;
import java.util.*;
/** Expression-trait association scan; this is not genetically predicted-expression TWAS. */
public class Twas {
  public static class Result {public int gene,n;public double beta,rSquared,pValue;}
  public static List<Result> scan(double[][] expression,double[] phenotype){
    if(expression.length!=phenotype.length)throw new IllegalArgumentException("sample rows mismatch");
    List<Result> out=new ArrayList<>();int n=phenotype.length;
    for(int g=0;g<expression[0].length;g++){
      double sx=0,sy=0,sxx=0,syy=0,sxy=0;int c=0;
      for(int i=0;i<n;i++){double x=expression[i][g],y=phenotype[i];if(Double.isFinite(x)&&Double.isFinite(y)){sx+=x;sy+=y;sxx+=x*x;syy+=y*y;sxy+=x*y;c++;}}
      Result r=new Result();r.gene=g;r.n=c;double dx=c*sxx-sx*sx,dy=c*syy-sy*sy;
      if(c<5||dx<=0||dy<=0){r.beta=r.rSquared=r.pValue=Double.NaN;}
      else{r.beta=(c*sxy-sx*sy)/dx;double corr=(c*sxy-sx*sy)/Math.sqrt(dx*dy);r.rSquared=corr*corr;
        double t=Math.abs(corr)*Math.sqrt((c-2)/Math.max(1e-15,1-corr*corr));r.pValue=2*(1-normalCdf(t));}
      out.add(r);
    }return out;
  }
  private static double normalCdf(double x){double z=x/Math.sqrt(2),t=1/(1+0.3275911*z);
    double erf=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-z*z);return .5*(1+erf);}
}
