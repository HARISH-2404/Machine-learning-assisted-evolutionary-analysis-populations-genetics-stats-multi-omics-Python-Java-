package org.ipgml.quantitativegenetics;
import java.util.*;
public final class CombiningAbility {
 public record Cross(String parent1,String parent2,double mean){}
 public static Map<String,Object> estimate(List<Cross> crosses){
  if(crosses.isEmpty())throw new IllegalArgumentException("No crosses");Set<String> ps=new TreeSet<>();double grand=0;for(Cross c:crosses){ps.add(c.parent1());ps.add(c.parent2());grand+=c.mean();}grand/=crosses.size();
  Map<String,Double> gca=new TreeMap<>();for(String p:ps){double sum=0;int n=0;for(Cross c:crosses)if(c.parent1().equals(p)||c.parent2().equals(p)){sum+=c.mean();n++;}gca.put(p,n==0?0:sum/n-grand);}
  List<Double> sca=new ArrayList<>();for(Cross c:crosses)sca.add(c.mean()-grand-gca.get(c.parent1())-gca.get(c.parent2()));
  Map<String,Object> out=new LinkedHashMap<>();out.put("grand_mean",grand);out.put("gca",gca);out.put("sca",sca);return out;
 }
}
