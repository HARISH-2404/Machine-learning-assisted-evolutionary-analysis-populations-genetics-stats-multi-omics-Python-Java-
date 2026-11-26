package org.ipgml.genomics;
import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/** Shared FASTA/VCF readers. Coordinates exposed by this package are 1-based inclusive. */
public final class GenomicsIO {
  private GenomicsIO() {}
  public static BufferedReader reader(String path) throws IOException {
    InputStream in=new FileInputStream(path);
    if(path.endsWith(".gz")) in=new GZIPInputStream(in);
    return new BufferedReader(new InputStreamReader(in));
  }
  public static LinkedHashMap<String,String> readFasta(String path) throws IOException {
    LinkedHashMap<String,String> m=new LinkedHashMap<>(); String id=null; StringBuilder seq=new StringBuilder();
    try(BufferedReader br=reader(path)){ String line;
      while((line=br.readLine())!=null){ line=line.trim(); if(line.isEmpty()) continue;
        if(line.startsWith(">")){ if(id!=null)m.put(id,seq.toString().toUpperCase(Locale.ROOT));
          String h=line.substring(1).trim(); if(h.isEmpty())throw new IOException("Empty FASTA header");
          id=h.split("\\s+")[0]; seq=new StringBuilder();
        } else { if(id==null)throw new IOException("Sequence before FASTA header"); seq.append(line); }
      }
    }
    if(id!=null)m.put(id,seq.toString().toUpperCase(Locale.ROOT)); return m;
  }
  public static int dosage(String gt) {
    if(gt==null||gt.isEmpty()||gt.contains("."))return -1;
    String[] a=gt.replace('|','/').split("/");
    if(a.length!=2)return -1;
    try{return Integer.parseInt(a[0])+Integer.parseInt(a[1]);}catch(NumberFormatException e){return -1;}
  }
  public static List<Map<String,Object>> readVcf(String path) throws IOException {
    List<Map<String,Object>> rows=new ArrayList<>();
    try(BufferedReader br=reader(path)){String line; List<String> samples=new ArrayList<>();
      while((line=br.readLine())!=null){
        if(line.startsWith("##"))continue;
        if(line.startsWith("#CHROM")){String[] h=line.split("\\t",-1);for(int i=9;i<h.length;i++)samples.add(h[i]);continue;}
        if(line.isEmpty()||line.startsWith("#"))continue;
        String[] f=line.split("\\t",-1);if(f.length<8)throw new IOException("Malformed VCF record");
        Map<String,Object> r=new LinkedHashMap<>();r.put("chrom",f[0]);r.put("pos",Integer.parseInt(f[1]));
        r.put("id",f[2]);r.put("ref",f[3]);r.put("alt",f[4]);r.put("qual",f[5]);r.put("filter",f[6]);
        List<Integer> ds=new ArrayList<>();
        if(f.length>8){String[] formats=f[8].split(":");int gi=Arrays.asList(formats).indexOf("GT");
          for(int i=9;i<f.length;i++){String[] vals=f[i].split(":",-1);ds.add(gi>=0&&gi<vals.length?dosage(vals[gi]):-1);}
        }
        r.put("samples",new ArrayList<>(samples));r.put("dosages",ds);rows.add(r);
      }
    }return rows;
  }
  public static int[] stats(String sequence) {
    int[] c=new int[6]; for(char x:sequence.toUpperCase(Locale.ROOT).toCharArray()){
      switch(x){case 'A':c[0]++;break;case 'C':c[1]++;break;case 'G':c[2]++;break;case 'T':c[3]++;break;case 'N':c[4]++;break;default:c[5]++;}
    }return c;
  }
}
