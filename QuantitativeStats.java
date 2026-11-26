package org.ipgml.quantitativegenetics;
import java.util.*;
public final class QuantitativeStats {
 public static Map<String,Double> summary(double[] x){if(x.length==0)throw new IllegalArgumentException("Empty values");double mean=Arrays.stream(x).average().orElseThrow(),ss=0;for(double v:x)ss+=(v-mean)*(v-mean);double sd=x.length>1?Math.sqrt(ss/(x.length-1)):Double.NaN;Map<String,Double>m=new LinkedHashMap<>();m.put("mean",mean);m.put("sd",sd);m.put("min",Arrays.stream(x).min().orElseThrow());m.put("max",Arrays.stream(x).max().orElseThrow());m.put("cv_percent",mean==0?Double.NaN:100*sd/mean);return m;}
 public static double broadSenseH2(double Vg,double Vp){if(Vg<0||Vp<=0)throw new IllegalArgumentException("Invalid variances");return Vg/Vp;}
 public static double narrowSenseH2(double Va,double Vp){if(Va<0||Vp<=0)throw new IllegalArgumentException("Invalid variances");return Va/Vp;}
 public static double geneticAdvance(double i,double h2,double sdP){if(i<0||h2<0||h2>1||sdP<0)throw new IllegalArgumentException("Invalid inputs");return i*h2*sdP;}
 public static double selectionResponse(double S,double h2){if(h2<0||h2>1)throw new IllegalArgumentException("h2 outside [0,1]");return S*h2;}
}
