package org.ipgml.moleculargenetics;
public final class MolecularGeneticsCli {
 public static void main(String[] args){if(args.length<2){System.err.println("Usage: reverse-complement|translate|gc|protein <sequence>");System.exit(2);}
 String cmd=args[0],seq=args[1];switch(cmd){case "reverse-complement"->System.out.println(SequenceUtils.reverseComplement(seq));case "translate"->System.out.println(Translation.translate(seq,0,false));case "gc"->System.out.println(SequenceUtils.gcPercent(seq));case "protein"->System.out.println(ProteinProperties.calculate(seq));default->throw new IllegalArgumentException("Unknown command "+cmd);}}
}
