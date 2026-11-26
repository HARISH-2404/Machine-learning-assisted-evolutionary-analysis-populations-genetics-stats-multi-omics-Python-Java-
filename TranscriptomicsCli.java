package org.ipgml.transcriptomics;
/** CLI for count-matrix QC and simple exploratory differential expression. */
public class TranscriptomicsCli {
  public static void main(String[] a)throws Exception{
    if(a.length<2){usage();return;}
    if(a[0].equals("qc")){
      if(a.length!=3)throw new IllegalArgumentException("qc counts.tsv qc.csv");
      ExpressionMatrix.readTsv(a[1]).writeQc(a[2]);
    } else if(a[0].equals("de")){
      if(a.length!=6)throw new IllegalArgumentException("de counts.tsv groups_csv control treatment output.csv");
      ExpressionMatrix m=ExpressionMatrix.readTsv(a[1]);String[] groups=a[2].split(",");
      DifferentialExpression.run(m,groups,a[3],a[4],a[5]);
    } else usage();
  }
  private static void usage(){System.out.println("TranscriptomicsCli qc counts.tsv qc.csv");System.out.println("TranscriptomicsCli de counts.tsv groups_csv control treatment output.csv");}
}
