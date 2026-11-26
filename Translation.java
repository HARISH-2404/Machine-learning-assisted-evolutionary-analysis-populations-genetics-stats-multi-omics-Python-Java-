package org.ipgml.moleculargenetics;
import java.util.*;
public final class Translation {
 private static final Map<String,Character> TABLE=new HashMap<>();
 static {
 String[] rows={"TTT F TTC F TTA L TTG L TCT S TCC S TCA S TCG S TAT Y TAC Y TAA * TAG * TGT C TGC C TGA * TGG W",
 "CTT L CTC L CTA L CTG L CCT P CCC P CCA P CCG P CAT H CAC H CAA Q CAG Q CGT R CGC R CGA R CGG R",
 "ATT I ATC I ATA I ATG M ACT T ACC T ACA T ACG T AAT N AAC N AAA K AAG K AGT S AGC S AGA R AGG R",
 "GTT V GTC V GTA V GTG V GCT A GCC A GCA A GCG A GAT D GAC D GAA E GAG E GGT G GGC G GGA G GGG G"};
 for(String row:rows){String[] t=row.split(" ");for(int i=0;i<t.length;i+=2)TABLE.put(t[i],t[i+1].charAt(0));}
 }
 public static String transcribe(String dna){return SequenceUtils.validateDna(dna).replace('T','U');}
 public static String translate(String sequence,int frame,boolean stopAtStop){
 if(frame<0||frame>2)throw new IllegalArgumentException("frame must be 0..2");
 String x=SequenceUtils.validateDna(sequence);StringBuilder p=new StringBuilder();
 for(int i=frame;i+2<x.length();i+=3){char aa=TABLE.getOrDefault(x.substring(i,i+3),'X');if(aa=='*'&&stopAtStop)break;p.append(aa);}return p.toString();}
}
