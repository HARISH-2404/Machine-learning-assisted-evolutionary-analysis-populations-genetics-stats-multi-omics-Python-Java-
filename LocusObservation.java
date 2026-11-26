package org.ipgml.integration;
/** Longitudinal allele-count observation at a genomic locus. */
public class LocusObservation {
 public final String chromosome,geneId,allele,population,location,timepoint;
 public final long position,generation,alleleCount,totalAlleles;
 public LocusObservation(String chromosome,long position,String geneId,String allele,String population,String location,String timepoint,long generation,long alleleCount,long totalAlleles){
  if(totalAlleles<=0 || alleleCount<0 || alleleCount>totalAlleles) throw new IllegalArgumentException("Invalid allele counts");
  this.chromosome=chromosome; this.position=position; this.geneId=geneId; this.allele=allele; this.population=population; this.location=location; this.timepoint=timepoint; this.generation=generation; this.alleleCount=alleleCount; this.totalAlleles=totalAlleles;
 }
 public double frequency(){return (double)alleleCount/totalAlleles;}
}
