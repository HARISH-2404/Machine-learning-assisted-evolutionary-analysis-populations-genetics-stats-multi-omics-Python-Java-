"""GFF3 gene models, coordinate validation and BED conversion."""
from .variant_annotation import read_gff_genes

def genes_to_bed(gff_path, bed_path):
    """Convert GFF3 gene rows (1-based inclusive) to BED (0-based half-open)."""
    genes=read_gff_genes(gff_path)
    with open(bed_path,"w",encoding="utf-8") as out:
        for g in genes:
            out.write(f'{g["chrom"]}\t{g["start"]-1}\t{g["end"]}\t{g["gene_id"]}\t0\t{g["strand"]}\n')
    return len(genes)

def genes_overlapping_interval(genes, chrom, start, end):
    """Return gene models overlapping a 1-based inclusive interval."""
    start=int(start); end=int(end)
    if start<1 or end<start: raise ValueError("invalid 1-based interval")
    return [g for g in genes if g["chrom"]==chrom and g["start"]<=end and g["end"]>=start]

def extract_gene_sequences(gff_path, fasta_path, output_fasta):
    """Extract gene-span sequences (not spliced transcripts); reverse-complement minus strand."""
    from ._core import read_fasta
    from .fasta import extract_interval
    genes=read_gff_genes(gff_path); genome=read_fasta(fasta_path)
    comp=str.maketrans("ACGTN","TGCAN")
    n=0
    with open(output_fasta,"w",encoding="utf-8") as out:
        for g in genes:
            if g["chrom"] not in genome: continue
            seq=genome[g["chrom"]][g["start"]-1:g["end"]]
            if g["strand"]=="-": seq=seq.translate(comp)[::-1]
            out.write(f'>{g["gene_id"]} {g["chrom"]}:{g["start"]}-{g["end"]}({g["strand"]})\n')
            for i in range(0,len(seq),60): out.write(seq[i:i+60]+"\n")
            n+=1
    return n
