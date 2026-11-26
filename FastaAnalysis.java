package org.ipgml.genomics;
import java.io.*;
import java.util.*;

/** FASTA sequence statistics and 1-based inclusive sequence extraction. */
public class FastaAnalysis {
  public static Map<String,Double> statistics(String sequence){
    int[] c=GenomicsIO.stats(sequence);int valid=c[0]+c[1]+c[2]+c[3];
    Map<String,Double> r=new LinkedHashMap<>();r.put("length",(double)sequence.length());
    r.put("A",(double)c[0]);r.put("C",(double)c[1]);r.put("G",(double)c[2]);r.put("T",(double)c[3]);r.put("N",(double)c[4]);
    r.put("GC_percent",valid==0?Double.NaN:100.0*(c[1]+c[2])/valid);return r;
  }
  public static String extract(String path,String id,int start,int end)throws IOException{
    if(start<1||end<start)throw new IllegalArgumentException("Require 1-based inclusive interval");
    String s=GenomicsIO.readFasta(path).get(id);if(s==null)throw new IllegalArgumentException("Sequence not found: "+id);
    if(end>s.length())throw new IllegalArgumentException("End exceeds sequence length");
    return s.substring(start-1,end);
  }
  public static void writeStatistics(String fasta,String output)throws IOException{
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      w.write("sequence_id,length,A,C,G,T,N,GC_percent\n");
      for(Map.Entry<String,String> e:GenomicsIO.readFasta(fasta).entrySet()){
        Map<String,Double>s=statistics(e.getValue());
        w.write(e.getKey()+","+s.get("length")+","+s.get("A")+","+s.get("C")+","+s.get("G")+","+s.get("T")+","+s.get("N")+","+s.get("GC_percent")+"\n");
      }
    }
  }
}
