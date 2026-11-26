package org.ipgml.populationgenetics;

import java.io.*;
import java.nio.file.*;
import java.util.*;

/** In-memory diploid biallelic genotype matrix. Dosage: 0,1,2; missing: NaN. */
public final class GenotypeData {
    public final List<String> sampleIds, populations, markers;
    public final double[][] dosage;

    public GenotypeData(List<String> sampleIds, List<String> populations,
                        List<String> markers, double[][] dosage) {
        this.sampleIds=sampleIds; this.populations=populations;
        this.markers=markers; this.dosage=dosage;
        if (dosage.length!=sampleIds.size() || populations.size()!=sampleIds.size())
            throw new IllegalArgumentException("Sample metadata dimensions do not match genotype rows.");
        for (double[] row:dosage) {
            if (row.length!=markers.size()) throw new IllegalArgumentException("Marker dimension mismatch.");
            for(double v:row) if(!Double.isNaN(v) && v!=0 && v!=1 && v!=2)
                throw new IllegalArgumentException("Dosage must be 0,1,2 or NaN.");
        }
    }

    /** Reads simple CSV: genotype_id,population_id,SNP_001,... (no quoted commas). */
    public static GenotypeData readCsv(Path path) throws IOException {
        List<String> lines=Files.readAllLines(path);
        if(lines.size()<2) throw new IllegalArgumentException("CSV must contain header and data.");
        String[] h=lines.get(0).split(",",-1);
        int id=-1,pop=-1; List<Integer> mi=new ArrayList<>(); List<String> markers=new ArrayList<>();
        for(int j=0;j<h.length;j++){
            if(h[j].equals("genotype_id")) id=j;
            else if(h[j].equals("population_id")) pop=j;
            else {mi.add(j); markers.add(h[j]);}
        }
        if(id<0 || pop<0 || markers.isEmpty())
            throw new IllegalArgumentException("CSV needs genotype_id, population_id and marker columns.");
        List<String> ids=new ArrayList<>(), pops=new ArrayList<>();
        double[][] x=new double[lines.size()-1][markers.size()];
        Set<String> seen=new HashSet<>();
        for(int i=1;i<lines.size();i++){
            String[] a=lines.get(i).split(",",-1);
            if(a.length!=h.length) throw new IllegalArgumentException("Malformed CSV row "+(i+1));
            if(!seen.add(a[id])) throw new IllegalArgumentException("Duplicate genotype_id: "+a[id]);
            ids.add(a[id]); pops.add(a[pop]);
            for(int k=0;k<mi.size();k++){
                String v=a[mi.get(k)].trim();
                x[i-1][k]=v.isEmpty()||v.equalsIgnoreCase("NA")||v.equalsIgnoreCase("NaN")
                    ? Double.NaN : Double.parseDouble(v);
            }
        }
        return new GenotypeData(ids,pops,markers,x);
    }

    public List<String> uniquePopulations(){
        return new ArrayList<>(new TreeSet<>(populations));
    }
    public static void writeCsv(Path path, String header, List<String> rows) throws IOException {
        Files.createDirectories(path.toAbsolutePath().getParent());
        try(BufferedWriter w=Files.newBufferedWriter(path)){
            w.write(header); w.newLine();
            for(String r:rows){w.write(r);w.newLine();}
        }
    }
    public static String fmt(double x){return Double.isFinite(x)?Double.toString(x):"";}
}
