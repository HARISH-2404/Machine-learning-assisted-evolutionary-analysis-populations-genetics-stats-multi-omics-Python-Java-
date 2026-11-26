package org.ipgml.integration;
import java.util.*;
/** Sample-by-feature numeric block. */
public class OmicsMatrix {
 public final List<String> samples, features; public final double[][] values;
 public OmicsMatrix(List<String> samples,List<String> features,double[][] values){
  if(values.length!=samples.size()) throw new IllegalArgumentException("Rows must match samples");
  for(double[] r:values) if(r.length!=features.size()) throw new IllegalArgumentException("Columns must match features");
  if(new HashSet<>(samples).size()!=samples.size() || new HashSet<>(features).size()!=features.size()) throw new IllegalArgumentException("IDs must be unique");
  this.samples=List.copyOf(samples); this.features=List.copyOf(features); this.values=values;
 }
}
