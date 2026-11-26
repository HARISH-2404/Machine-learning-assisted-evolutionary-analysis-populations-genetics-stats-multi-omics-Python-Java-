package org.ipgml.moleculargenetics;
import java.util.*;
public final class SequenceUtils {
 private SequenceUtils(){}
 public static String validateDna(String s){String x=s.replaceAll("\\s+","").toUpperCase().replace('U','T');
 if(x.isEmpty() || !x.matches("[ACGTN]+")) throw new IllegalArgumentException("Invalid/empty DNA sequence"); return x;}
 public static String reverseComplement(String s){String x=validateDna(s); StringBuilder b=new StringBuilder();
 for(int i=x.length()-1;i>=0;i--) b.append(switch(x.charAt(i)){case 'A'->'T';case 'T'->'A';case 'C'->'G';case 'G'->'C';default->'N';}); return b.toString();}
 public static double gcPercent(String s){String x=validateDna(s);int gc=0,n=0;for(char c:x.toCharArray()){if(c!='N'){n++;if(c=='G'||c=='C')gc++;}}return n==0?Double.NaN:100.0*gc/n;}
}
