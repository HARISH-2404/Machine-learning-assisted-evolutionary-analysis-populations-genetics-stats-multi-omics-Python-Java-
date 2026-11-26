package org.ipgml.moleculargenetics;
import java.nio.file.*;import java.io.*;import java.util.*;
/** Reads a simple TSV: gene_id, description, go_terms, kegg_pathways. */
public final class FunctionalAnnotation {
 public static Map<String,Map<String,String>> loadTsv(Path path)throws IOException{
  List<String> lines=Files.readAllLines(path);if(lines.isEmpty())throw new IllegalArgumentException("Empty annotation file");
  String[] h=lines.get(0).split("\t",-1);Map<String,Integer> ix=new HashMap<>();for(int i=0;i<h.length;i++)ix.put(h[i].trim(),i);
  if(!ix.containsKey("gene_id"))throw new IllegalArgumentException("gene_id column required");
  Map<String,Map<String,String>> out=new LinkedHashMap<>();
  for(int r=1;r<lines.size();r++){String[] a=lines.get(r).split("\t",-1);String id=a[ix.get("gene_id")].trim();if(id.isEmpty())continue;
   Map<String,String> row=new LinkedHashMap<>();for(String col:List.of("description","go_terms","kegg_pathways"))row.put(col,ix.containsKey(col)&&ix.get(col)<a.length?a[ix.get(col)].trim():"");out.put(id,row);}
  return out;
 }
 public static double coverage(Collection<String> genes,Map<String,?> annotations){if(genes.isEmpty())return 0;long n=genes.stream().distinct().filter(annotations::containsKey).count();return 100.0*n/genes.stream().distinct().count();}
}
