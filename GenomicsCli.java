package org.ipgml.genomics;
/** Command line entry point. Compile all sources, then run one of the commands below. */
public class GenomicsCli {
  public static void main(String[] a)throws Exception{
    if(a.length<2){usage();return;}
    switch(a[0]){
      case "fasta-stats":
        if(a.length!=3)throw new IllegalArgumentException("fasta-stats input.fa output.csv");
        FastaAnalysis.writeStatistics(a[1],a[2]);break;
      case "vcf-stats":
        if(a.length!=3)throw new IllegalArgumentException("vcf-stats input.vcf[.gz] output.csv");
        VariantAnalysis.summarizeVcf(a[1],a[2]);break;
      case "variant-density":
        if(a.length!=5)throw new IllegalArgumentException("variant-density input.vcf window_size step output.csv");
        GenomeWindows.variantDensity(a[1],Integer.parseInt(a[2]),Integer.parseInt(a[3]),a[4]);break;
      case "annotate":
        if(a.length!=5)throw new IllegalArgumentException("annotate input.vcf genes.gff3 upstream_bp output.csv");
        GffAnnotation.annotateVcf(a[1],a[2],a[4],Integer.parseInt(a[3]));break;
      case "genes-bed":
        if(a.length!=3)throw new IllegalArgumentException("genes-bed genes.gff3 output.bed");
        System.out.println("Genes converted: "+GeneCoordinates.toBed(a[1],a[2]));break;
      case "gene-sequences":
        if(a.length!=4)throw new IllegalArgumentException("gene-sequences genes.gff3 genome.fa output.fa");
        System.out.println("Gene spans extracted: "+GeneCoordinates.extractGeneSpans(a[1],a[2],a[3]));break;
      default:usage();
    }
  }
  private static void usage(){
    System.out.println("GenomicsCli commands:");
    System.out.println(" fasta-stats input.fa output.csv");
    System.out.println(" vcf-stats input.vcf[.gz] output.csv");
    System.out.println(" variant-density input.vcf window_size step output.csv");
    System.out.println(" annotate input.vcf genes.gff3 upstream_bp output.csv");
    System.out.println(" genes-bed genes.gff3 output.bed");
    System.out.println(" gene-sequences genes.gff3 genome.fa output.fa");
  }
}
