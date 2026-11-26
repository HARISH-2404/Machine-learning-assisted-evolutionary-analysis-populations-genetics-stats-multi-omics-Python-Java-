package org.ipgml.quantitativegenetics;
import java.util.*;
public final class GxeStability {
 public static Map<String,double[]> finlayWilkinson(double[][] genotypeByEnvironment){
  int g=genotypeByEnvironment.length;if(g==0)throw new IllegalArgumentException("No genotypes");int e=genotypeByEnvironment[0].length;if(e<3)throw new IllegalArgumentException("At least 3 environments");
  double[] env=new double[e];for(int j=0;j<e;j++){for(int i=0;i<g;i++)env[j]+=genotypeByEnvironment[i][j];env[j]/=g;}
  double em=Arrays.stream(env).average().orElseThrow();for(int j=0;j<e;j++)env[j]-=em;
  double sx=0;for(double x:env)sx+=x*x;if(sx==0)throw new IllegalArgumentException("No environment variation");
  Map<String,double[]> out=new LinkedHashMap<>();for(int i=0;i<g;i++){double ym=Arrays.stream(genotypeByEnvironment[i]).average().orElseThrow(),num=0,den=0;for(int j=0;j<e;j++){num+=env[j]*(genotypeByEnvironment[i][j]-ym);den+=env[j]*env[j];}double slope=num/den,intercept=ym;double sse=0,sst=0;for(int j=0;j<e;j++){double d=genotypeByEnvironment[i][j]-ym,p=ym+slope*env[j];sse+=(genotypeByEnvironment[i][j]-p)*(genotypeByEnvironment[i][j]-p);sst+=d*d;}out.put("G"+(i+1),new double[]{intercept,slope,sst==0?Double.NaN:1-sse/sst});}return out;
 }
}
