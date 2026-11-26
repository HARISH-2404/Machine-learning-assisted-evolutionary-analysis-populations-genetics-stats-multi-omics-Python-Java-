package org.ipgml.moleculargenetics;
import java.util.*;
public final class MotifScanner {
 private static final Map<Character,String> I=new HashMap<>();
 static {I.put('A',"A");I.put('C',"C");I.put('G',"G");I.put('T',"T");I.put('R',"[AG]");I.put('Y',"[CT]");I.put('S',"[GC]");I.put('W',"[AT]");I.put('K',"[GT]");I.put('M',"[AC]");I.put('B',"[CGT]");I.put('D',"[AGT]");I.put('H',"[ACT]");I.put('V',"[ACG]");I.put('N',"[ACGTN]");}
 public static List<Integer> scan(String sequence,String motif){String s=SequenceUtils.validateDna(sequence),m=motif.toUpperCase().replace('U','T');StringBuilder r=new StringBuilder();for(char c:m.toCharArray()){if(!I.containsKey(c))throw new IllegalArgumentException("Invalid IUPAC motif");r.append(I.get(c));}
 List<Integer> hits=new ArrayList<>();for(int i=0;i<=s.length()-m.length();i++)if(s.substring(i,i+m.length()).matches(r.toString()))hits.add(i+1);return hits;}
}
