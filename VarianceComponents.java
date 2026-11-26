package org.ipgml.quantitativegenetics;
import java.util.*;
public final class VarianceComponents {
 public static Map<String,Double> oneWay(Map<String,double[]> groups){
  if(groups.size()<2)throw new IllegalArgumentException("At least two groups");int r=groups.values().iterator().next().length;if(r<2)throw new IllegalArgumentException("At least two replicates");
  double grand=0;int count=0;Map<String,Double> means=new LinkedHashMap<>();
  for(var e:groups.entrySet()){if(e.getValue().length!=r)throw new IllegalArgumentException("Balanced replication required");double m=Arrays.stream(e.getValue()).average().orElseThrow();means.put(e.getKey(),m);grand+=Arrays.stream(e.getValue()).sum();count+=r;}
  grand/=count;double ssb=0,ssw=0;for(var e:groups.entrySet()){double m=means.get(e.getKey());ssb+=r*(m-grand)*(m-grand);for(double v:e.getValue())ssw+=(v-m)*(v-m);}
  double msb=ssb/(groups.size()-1),msw=ssw/(groups.size()*(r-1));Map<String,Double> out=new LinkedHashMap<>();out.put("genetic_variance",Math.max(0,(msb-msw)/r));out.put("residual_variance",msw);out.put("ms_genotype",msb);out.put("ms_error",msw);return out;
 }
}
