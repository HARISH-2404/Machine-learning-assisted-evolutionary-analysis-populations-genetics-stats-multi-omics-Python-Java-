package org.ipgml.forecasting;
public class ForecastCli {
 public static void main(String[] args){
  double[][] paths=WrightFisher.simulate(0.15,20,200,0.05,0.5,1000,2026L);
  double[] s=ForecastSummary.meanAndSd(paths,20);
  System.out.printf("Generation 20 mean frequency=%.4f SD=%.4f P(freq>=0.8)=%.4f%n",
    s[0],s[1],ForecastSummary.thresholdProbability(paths,20,0.8));
 }
}
