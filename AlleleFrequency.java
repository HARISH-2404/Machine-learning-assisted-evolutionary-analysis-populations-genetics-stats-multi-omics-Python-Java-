package org.ipgml.populationgenetics;
import java.util.*;
import java.nio.file.*;
/** Estimates alternate allele frequency p, reference frequency q, MAF and call rate. */
public class AlleleFrequency {
    public static List<String> estimate(GenotypeData d) {
        List<String> out=new ArrayList<>();
        for(String pop:d.uniquePopulations()) for(int m=0;m<d.markers.size();m++){
            double sum=0; int n=0, total=0;
            for(int i=0;i<d.sampleIds.size();i++) if(d.populations.get(i).equals(pop)){
                total++; double v=d.dosage[i][m]; if(Double.isFinite(v)){sum+=v;n++;}
            }
            double p=n>0?sum/(2*n):Double.NaN, q=1-p;
            out.add(pop+","+d.markers.get(m)+","+n+","+GenotypeData.fmt(total>0?(double)n/total:Double.NaN)+","+
                GenotypeData.fmt(p)+","+GenotypeData.fmt(q)+","+GenotypeData.fmt(Math.min(p,q)));
        }
        return out;
    }
    public static void write(GenotypeData d,Path out)throws Exception{
        GenotypeData.writeCsv(out,"population_id,marker_id,n_called,call_rate,alt_frequency,ref_frequency,maf",estimate(d));
    }
}
