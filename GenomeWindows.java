package org.ipgml.genomics;
import java.io.*;
import java.util.*;

/** Sliding-window genomic summaries. Coordinates are 1-based inclusive. */
public class GenomeWindows {
  public static List<int[]> windows(int length,int windowSize,int step){
    if(length<1||windowSize<1||step<1)throw new IllegalArgumentException("positive lengths required");
    List<int[]> out=new ArrayList<>();
    for(int start=1;start<=length;start+=step){int end=Math.min(length,start+windowSize-1);out.add(new int[]{start,end});if(end==length)break;}
    return out;
  }
  public static void variantDensity(String vcf,int windowSize,int step,String output)throws IOException{
    Map<String,List<Integer>> byChr=new LinkedHashMap<>();
    for(Map<String,Object> r:GenomicsIO.readVcf(vcf)){
      String chr=(String)r.get("chrom");byChr.computeIfAbsent(chr,k->new ArrayList<>()).add((Integer)r.get("pos"));
    }
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      w.write("chrom,start,end,variant_count,density_per_kb\n");
      for(String chr:byChr.keySet()){
        List<Integer> p=byChr.get(chr);int max=Collections.max(p);
        for(int[] win:windows(max,windowSize,step)){int n=0;for(int x:p)if(x>=win[0]&&x<=win[1])n++;
          w.write(chr+","+win[0]+","+win[1]+","+n+","+(n*1000.0/(win[1]-win[0]+1))+"\n");
        }
      }
    }
  }
}
