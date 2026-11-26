package org.ipgml.moleculargenetics;
import java.util.*;
/** Directed regulator-to-target graph utilities. */
public final class RegulatoryNetwork {
 public record Edge(String regulator,String target,String evidence,double score){}
 public static Map<String,Object> build(List<Edge> edges,double minimumScore){
  Set<String> nodes=new TreeSet<>();List<Edge> kept=new ArrayList<>();for(Edge e:edges)if(e.score()>=minimumScore){nodes.add(e.regulator());nodes.add(e.target());kept.add(e);}
  Map<String,Object> out=new LinkedHashMap<>();out.put("nodes",nodes);out.put("edges",kept);return out;
 }
 public static Map<String,Integer> outDegree(List<Edge> edges){Map<String,Integer>d=new TreeMap<>();for(Edge e:edges)d.merge(e.regulator(),1,Integer::sum);return d;}
 public static Set<String> reachable(List<Edge> edges,String start){Map<String,List<String>> a=new HashMap<>();for(Edge e:edges)a.computeIfAbsent(e.regulator(),k->new ArrayList<>()).add(e.target());Set<String> seen=new TreeSet<>();Deque<String> q=new ArrayDeque<>(a.getOrDefault(start,List.of()));while(!q.isEmpty()){String x=q.remove();if(seen.add(x))q.addAll(a.getOrDefault(x,List.of()));}return seen;}
}
