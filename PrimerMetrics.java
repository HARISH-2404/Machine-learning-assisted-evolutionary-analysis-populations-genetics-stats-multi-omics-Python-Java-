package org.ipgml.moleculargenetics;
import java.util.*;
public final class PrimerMetrics {
 public static Map<String,Object> calculate(String primer){String s=SequenceUtils.validateDna(primer);if(s.contains("N"))throw new IllegalArgumentException("Primer cannot contain N");
 int at=0,gc=0,maxRun=1,run=1;for(int i=0;i<s.length();i++){char c=s.charAt(i);if(c=='A'||c=='T')at++;else gc++;if(i>0){if(c==s.charAt(i-1)){run++;maxRun=Math.max(maxRun,run);}else run=1;}}
 Map<String,Object> m=new LinkedHashMap<>();m.put("sequence",s);m.put("length",s.length());m.put("gc_percent",100.0*gc/s.length());m.put("tm_wallace_c",2*at+4*gc);m.put("reverse_complement",SequenceUtils.reverseComplement(s));m.put("homopolymer_max",maxRun);return m;}
}
