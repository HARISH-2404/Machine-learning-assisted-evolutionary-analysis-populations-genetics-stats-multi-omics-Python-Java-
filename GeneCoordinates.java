package org.ipgml.genomics;
import java.io.*;
import java.util.*;

/** GFF3 gene-span extraction and BED coordinate conversion. */
public class GeneCoordinates {
  public static int toBed(String gff,String bed)throws IOException{
    List<GffAnnotation.Gene> genes=GffAnnotation.readGenes(gff);
    try(BufferedWriter w=new BufferedWriter(new FileWriter(bed))){
      for(GffAnnotation.Gene g:genes)w.write(g.chrom+"\t"+(g.start-1)+"\t"+g.end+"\t"+g.id+"\t0\t"+g.strand+"\n");
    }return genes.size();
  }
  public static int extractGeneSpans(String gff,String fasta,String output)throws IOException{
    List<GffAnnotation.Gene> genes=GffAnnotation.readGenes(gff);
    Map<String,String> genome=GenomicsIO.readFasta(fasta);int n=0;
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      for(GffAnnotation.Gene g:genes){String chr=genome.get(g.chrom);if(chr==null||g.end>chr.length())continue;
        String s=chr.substring(g.start-1,g.end);
        if(g.strand.equals("-"))s=reverseComplement(s);
        w.write(">"+g.id+" "+g.chrom+":"+g.start+"-"+g.end+"("+g.strand+")\n");
        for(int i=0;i<s.length();i+=60)w.write(s.substring(i,Math.min(i+60,s.length()))+"\n");n++;
      }
    }return n;
  }
  private static String reverseComplement(String s){
    StringBuilder b=new StringBuilder();for(int i=s.length()-1;i>=0;i--){char c=s.charAt(i);
      switch(c){case 'A':b.append('T');break;case 'T':b.append('A');break;case 'C':b.append('G');break;case 'G':b.append('C');break;default:b.append('N');}
    }return b.toString();
  }
}
