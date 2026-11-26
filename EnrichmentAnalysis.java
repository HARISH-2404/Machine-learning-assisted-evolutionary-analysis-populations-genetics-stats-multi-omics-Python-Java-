package org.ipgml.moleculargenetics;
import java.util.*;
/** Hypergeometric over-representation with BH adjustment; mapping is term -> member genes. */
public final class EnrichmentAnalysis {
 private static double logChoose(int n,int k){if(k<0||k>n)return Double.NEGATIVE_INFINITY;k=Math.min(k,n-k);double x=0;for(int i=1;i<=k;i++)x+=Math.log(n-k+i)-Math.log(i);return x;}
 private static double tail(int M,int n,int N,int k){double sum=0;int hi=Math.min(n,N),lo=Math.max(k,N-(M-n));for(int x=lo;x<=hi;x++)sum+=Math.exp(logChoose(n,x)+logChoose(M-n,N-x)-logChoose(M,N));return Math.min(1,sum);}
 public static List<Map<String,Object>> overRepresentation(Collection<String> selected,Map<String,? extends Collection<String>> termGenes,Collection<String> background,int minSize){
  Set<String> bg=new HashSet<>(background),sel=new HashSet<>(selected);sel.retainAll(bg);int M=bg.size(),N=sel.size();List<Map<String,Object>> out=new ArrayList<>();if(M==0||N==0)return out;
  for(var entry:termGenes.entrySet()){Set<String> members=new HashSet<>(entry.getValue());members.retainAll(bg);int n=members.size();if(n<minSize)continue;Set<String> overlap=new TreeSet<>(members);overlap.retainAll(sel);int k=overlap.size();if(k==0)continue;
   Map<String,Object> row=new LinkedHashMap<>();row.put("term",entry.getKey());row.put("overlap",k);row.put("term_size",n);row.put("genes",overlap);row.put("fold_enrichment",(k/(double)N)/(n/(double)M));row.put("p_value",tail(M,n,N,k));out.add(row);}
  out.sort(Comparator.comparingDouble(r->(double)r.get("p_value")));int m=out.size();double running=1;for(int i=m-1;i>=0;i--){double p=(double)out.get(i).get("p_value");running=Math.min(running,p*m/(i+1));out.get(i).put("p_adjusted",Math.min(1,running));}return out;
 }
}
