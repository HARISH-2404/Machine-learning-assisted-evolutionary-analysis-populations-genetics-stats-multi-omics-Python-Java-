package org.ipgml.genomics;
import java.io.*;
import java.util.*;

/** Basic GFF3 gene-span annotation; does not predict coding consequences. */
public class GffAnnotation {
  public static class Gene {
    public String chrom,id,strand; public int start,end;
    Gene(String c,int s,int e,String id,String strand){chrom=c;start=s;end=e;this.id=id;this.strand=strand;}
  }
  public static List<Gene> readGenes(String gff)throws IOException{
    List<Gene> genes=new ArrayList<>();
    try(BufferedReader br=new BufferedReader(new FileReader(gff))){String line;
      while((line=br.readLine())!=null){if(line.isEmpty()||line.startsWith("#"))continue;String[] f=line.split("\\t");
        if(f.length!=9||!f[2].equalsIgnoreCase("gene"))continue;
        String id="unknown";for(String a:f[8].split(";"))if(a.startsWith("ID=")){id=a.substring(3);break;}
        genes.add(new Gene(f[0],Integer.parseInt(f[3]),Integer.parseInt(f[4]),id,f[6]));
      }
    }return genes;
  }
  public static String[] annotate(String chrom,int pos,List<Gene> genes,int upstream){
    Gene best=null;int bestDist=Integer.MAX_VALUE;String region="intergenic";
    for(Gene g:genes)if(g.chrom.equals(chrom)){
      if(pos>=g.start&&pos<=g.end){int d=0;if(d<bestDist){best=g;bestDist=d;region="genic";}}
      else {int d=Math.min(Math.abs(pos-g.start),Math.abs(pos-g.end));if(d<=upstream&&d<bestDist){best=g;bestDist=d;region="near_gene";}}
    }
    return best==null?new String[]{"","intergenic",""}:new String[]{best.id,region,Integer.toString(bestDist)};
  }
  public static void annotateVcf(String vcf,String gff,String output,int upstream)throws IOException{
    List<Gene> genes=readGenes(gff);
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      w.write("chrom,pos,id,gene_id,region,distance\n");
      for(Map<String,Object> v:GenomicsIO.readVcf(vcf)){
        String[] a=annotate((String)v.get("chrom"),(Integer)v.get("pos"),genes,upstream);
        w.write(v.get("chrom")+","+v.get("pos")+","+v.get("id")+","+a[0]+","+a[1]+","+a[2]+"\n");
      }
    }
  }
}
