package org.ipgml.moleculargenetics;
import java.util.*;
public final class GeneStructure {
 public static Map<String,Object> summarize(List<int[]> exons){if(exons==null||exons.isEmpty())throw new IllegalArgumentException("Exons required");
 List<int[]> x=new ArrayList<>();for(int[] e:exons){if(e.length!=2||e[0]<1||e[1]<e[0])throw new IllegalArgumentException("Invalid exon");x.add(e.clone());}
 x.sort(Comparator.comparingInt(a->a[0]));int spliced=0;List<Integer> lens=new ArrayList<>();List<String> introns=new ArrayList<>();
 for(int i=0;i<x.size();i++){int len=x.get(i)[1]-x.get(i)[0]+1;spliced+=len;lens.add(len);if(i>0){if(x.get(i)[0]<=x.get(i-1)[1])throw new IllegalArgumentException("Overlapping exons");introns.add((x.get(i-1)[1]+1)+"-"+(x.get(i)[0]-1));}}
 Map<String,Object> out=new LinkedHashMap<>();out.put("exon_count",x.size());out.put("exon_lengths",lens);out.put("intron_intervals",introns);out.put("spliced_length",spliced);out.put("genomic_span",x.get(x.size()-1)[1]-x.get(0)[0]+1);return out;}
}
