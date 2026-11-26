package org.ipgml.transcriptomics;
import java.io.*;
import java.util.*;

/** GFF3 transcript-to-gene mapping and feature parsing. */
public class TranscriptAnnotation {
  public static Map<String,String> transcriptGeneMap(String gff)throws IOException{
    Map<String,String> map=new LinkedHashMap<>();
    try(BufferedReader br=new BufferedReader(new FileReader(gff))){String line;
      while((line=br.readLine())!=null){if(line.isEmpty()||line.startsWith("#"))continue;String[] f=line.split("\\t");
        if(f.length!=9||!(f[2].equalsIgnoreCase("mRNA")||f[2].equalsIgnoreCase("transcript")))continue;
        String id=null,parent=null;for(String a:f[8].split(";")){if(a.startsWith("ID="))id=a.substring(3);if(a.startsWith("Parent="))parent=a.substring(7);}
        if(id!=null&&parent!=null)map.put(id,parent);
      }
    }return map;
  }
}
