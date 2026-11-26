package org.ipgml.genomics;
import java.io.*;
import java.util.*;

/** Diploid biallelic VCF dosage summaries (0/1/2); missing genotypes are -1. */
public class VariantAnalysis {
  public static Map<String,Double> summarize(List<Integer> dosage){
    int called=0,het=0;double alt=0;
    for(int d:dosage)if(d>=0&&d<=2){called++;alt+=d;if(d==1)het++;}
    double p=called==0?Double.NaN:alt/(2.0*called);
    Map<String,Double> r=new LinkedHashMap<>();
    r.put("n_called",(double)called);r.put("n_missing",(double)(dosage.size()-called));
    r.put("call_rate",dosage.isEmpty()?Double.NaN:called/(double)dosage.size());
    r.put("alt_frequency",p);r.put("maf",Double.isNaN(p)?Double.NaN:Math.min(p,1-p));
    r.put("observed_heterozygosity",called==0?Double.NaN:het/(double)called);return r;
  }
  public static void summarizeVcf(String vcf,String output)throws IOException{
    try(BufferedWriter w=new BufferedWriter(new FileWriter(output))){
      w.write("chrom,pos,id,ref,alt,n_called,n_missing,call_rate,alt_frequency,maf,observed_heterozygosity\n");
      for(Map<String,Object> r:GenomicsIO.readVcf(vcf)){
        @SuppressWarnings("unchecked") List<Integer>d=(List<Integer>)r.get("dosages");
        Map<String,Double>s=summarize(d);
        w.write(r.get("chrom")+","+r.get("pos")+","+r.get("id")+","+r.get("ref")+","+r.get("alt")+","+s.get("n_called")+","+s.get("n_missing")+","+s.get("call_rate")+","+s.get("alt_frequency")+","+s.get("maf")+","+s.get("observed_heterozygosity")+"\n");
      }
    }
  }
}
