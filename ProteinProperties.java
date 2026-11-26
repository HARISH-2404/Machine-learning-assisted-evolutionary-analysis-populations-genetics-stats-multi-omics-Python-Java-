package org.ipgml.moleculargenetics;
import java.util.*;
public final class ProteinProperties {
 private static final Map<Character,Double> MASS=Map.ofEntries(Map.entry('A',89.09),Map.entry('R',174.20),Map.entry('N',132.12),Map.entry('D',133.10),Map.entry('C',121.15),Map.entry('E',147.13),Map.entry('Q',146.15),Map.entry('G',75.07),Map.entry('H',155.16),Map.entry('I',131.17),Map.entry('L',131.17),Map.entry('K',146.19),Map.entry('M',149.21),Map.entry('F',165.19),Map.entry('P',115.13),Map.entry('S',105.09),Map.entry('T',119.12),Map.entry('W',204.23),Map.entry('Y',181.19),Map.entry('V',117.15));
 public static Map<String,Object> calculate(String protein){String s=protein.replaceAll("\\s+","").toUpperCase();if(s.isEmpty())throw new IllegalArgumentException("Empty protein");double mass=0;Map<Character,Integer> counts=new TreeMap<>();for(char c:s.toCharArray()){if(!MASS.containsKey(c))throw new IllegalArgumentException("Invalid amino acid "+c);mass+=MASS.get(c);counts.merge(c,1,Integer::sum);}mass-=18.015*(s.length()-1);Map<String,Object> out=new LinkedHashMap<>();out.put("length_aa",s.length());out.put("molecular_mass_da_approx",mass);out.put("amino_acid_counts",counts);return out;}
}
