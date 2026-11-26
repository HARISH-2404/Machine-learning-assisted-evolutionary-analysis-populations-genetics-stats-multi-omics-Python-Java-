package org.ipgml.genomics;
import java.util.*;
/** Folded site-frequency spectrum and windowed expected heterozygosity. */
public class SelectionSweeps {
  public static Map<Integer,Integer> foldedSfs(double[][] dosage){
    Map<Integer,Integer> sfs=new TreeMap<>();
    for(int j=0;j<dosage[0].length;j++){int ac=0,an=0;
      for(double[] row:dosage)if(Double.isFinite(row[j])){ac+=(int)row[j];an+=2;}
      if(an>0){int minor=Math.min(ac,an-ac);sfs.put(minor,sfs.getOrDefault(minor,0)+1);}
    }return sfs;
  }
  public static class Window {public int start,end,nVariants;public double meanHe;}
  public static List<Window> heterozygosityWindows(double[][] dosage,int[] positions,int windowSize){
    if(dosage[0].length!=positions.length||windowSize<1)throw new IllegalArgumentException("dimensions/window size invalid");
    int min=Arrays.stream(positions).min().orElse(1),max=Arrays.stream(positions).max().orElse(1);
    List<Window> out=new ArrayList<>();
    for(int start=((min-1)/windowSize)*windowSize+1;start<=max;start+=windowSize){
      Window w=new Window();w.start=start;w.end=start+windowSize-1;double sum=0;
      for(int j=0;j<positions.length;j++)if(positions[j]>=w.start&&positions[j]<=w.end){
        int n=0,ac=0;for(double[] row:dosage)if(Double.isFinite(row[j])){n++;ac+=(int)row[j];}
        if(n>0){double p=ac/(2.0*n);sum+=2*p*(1-p);w.nVariants++;}
      }
      w.meanHe=w.nVariants>0?sum/w.nVariants:Double.NaN;out.add(w);
    }return out;
  }
}
