package org.ipgml.transcriptomics;
import java.util.*;
/** ASE imbalance using a two-sided normal approximation to binomial balance. */
public class AlleleSpecificExpression {
  public static class Result {public int ref,alt,total;public double altFraction,pValue;}
  public static List<Result> test(int[][] counts,int minTotal){
    List<Result> out=new ArrayList<>();
    for(int[] row:counts){if(row.length!=2||row[0]<0||row[1]<0)throw new IllegalArgumentException("rows must be [ref,alt] nonnegative counts");
      Result r=new Result();r.ref=row[0];r.alt=row[1];r.total=row[0]+row[1];r.altFraction=r.total==0?Double.NaN:row[1]/(double)r.total;
      if(r.total<minTotal)r.pValue=Double.NaN;
      else{double z=Math.abs(row[1]-r.total*.5)/Math.sqrt(r.total*.25);r.pValue=2*(1-normalCdf(z));}
      out.add(r);
    }return out;
  }
  private static double normalCdf(double x){double z=x/Math.sqrt(2),t=1/(1+0.3275911*z);
    double erf=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-z*z);return .5*(1+erf);}
}
