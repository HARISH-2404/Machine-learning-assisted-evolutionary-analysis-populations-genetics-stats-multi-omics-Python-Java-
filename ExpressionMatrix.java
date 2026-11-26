package org.ipgml.transcriptomics;
import java.io.*;
import java.util.*;

/** Dense genes-by-samples raw count matrix reader and validation. */
public class ExpressionMatrix {
  public final List<String> genes,samples;
  public final double[][] counts; // genes x samples
  private ExpressionMatrix(List<String> g,List<String>s,double[][]c){genes=g;samples=s;counts=c;}
  public static ExpressionMatrix readTsv(String path)throws IOException{
    try(BufferedReader br=new BufferedReader(new FileReader(path))){
      String h=br.readLine();if(h==null)throw new IOException("Empty count matrix");
      String[] head=h.split("\\t",-1);if(head.length<2)throw new IOException("Expected gene ID plus samples");
      List<String> samples=new ArrayList<>();for(int i=1;i<head.length;i++)samples.add(head[i]);
      List<String> genes=new ArrayList<>();List<double[]> rows=new ArrayList<>();Set<String> seen=new HashSet<>();String line;
      while((line=br.readLine())!=null){if(line.isEmpty())continue;String[] f=line.split("\\t",-1);
        if(f.length!=head.length)throw new IOException("Column count mismatch");
        if(!seen.add(f[0]))throw new IOException("Duplicate gene ID "+f[0]);
        double[] x=new double[samples.size()];
        for(int i=0;i<x.length;i++){x[i]=Double.parseDouble(f[i+1]);if(!Double.isFinite(x[i])||x[i]<0||x[i]!=Math.rint(x[i]))throw new IOException("Raw counts must be finite non-negative integers");}
        genes.add(f[0]);rows.add(x);
      }
      return new ExpressionMatrix(genes,samples,rows.toArray(new double[0][]));
    }
  }
  public double[] librarySizes(){double[] s=new double[samples.size()];for(double[] row:counts)for(int j=0;j<s.length;j++)s[j]+=row[j];return s;}
  public double[][] cpm(boolean log2,double prior){
    double[] libs=librarySizes();double[][] out=new double[counts.length][samples.size()];
    for(int i=0;i<counts.length;i++)for(int j=0;j<samples.size();j++){
      if(libs[j]<=0)throw new IllegalArgumentException("Zero-size library");
      double v=counts[i][j]*1e6/libs[j];out[i][j]=log2?Math.log(v+prior)/Math.log(2):v;
    }return out;
  }
  public void writeQc(String output)throws IOException{
    double[] lib=librarySizes();
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      w.write("sample_id,library_size,detected_genes,zero_fraction,median_gene_count\n");
      for(int j=0;j<samples.size();j++){int detected=0,zero=0;double[] vals=new double[counts.length];
        for(int i=0;i<counts.length;i++){vals[i]=counts[i][j];if(vals[i]>0)detected++;else zero++;}
        Arrays.sort(vals);double med=vals.length%2==1?vals[vals.length/2]:(vals[vals.length/2-1]+vals[vals.length/2])/2;
        w.write(samples.get(j)+","+lib[j]+","+detected+","+(zero/(double)counts.length)+","+med+"\n");
      }
    }
  }
}
