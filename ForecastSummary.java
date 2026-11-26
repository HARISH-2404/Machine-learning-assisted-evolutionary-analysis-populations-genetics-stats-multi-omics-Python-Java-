package org.ipgml.forecasting;
import java.util.*;
public final class ForecastSummary {
 private ForecastSummary(){}
 public static double[] meanAndSd(double[][] paths,int generation){
  double mean=0;for(double[] p:paths)mean+=p[generation]/paths.length;
  double v=0;for(double[] p:paths)v+=Math.pow(p[generation]-mean,2)/paths.length;
  return new double[]{mean,Math.sqrt(v)};
 }
 public static double thresholdProbability(double[][] paths,int generation,double threshold){
  int n=0;for(double[] p:paths)if(p[generation]>=threshold)n++;
  return (double)n/paths.length;
 }
}
