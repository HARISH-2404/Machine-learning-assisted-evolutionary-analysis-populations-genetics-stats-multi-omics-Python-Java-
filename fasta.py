"""FASTA sequence summaries and indexed interval extraction."""
from ._core import read_fasta, validate_interval

def sequence_statistics(sequence):
    s=sequence.upper(); valid=[x for x in s if x in "ACGT"]
    n=len(s); gc=sum(x in "GC" for x in valid)
    return {"length":n,"A":s.count("A"),"C":s.count("C"),"G":s.count("G"),
            "T":s.count("T"),"N":s.count("N"),
            "gc_percent":(100*gc/len(valid)) if valid else float("nan")}

def fasta_statistics(path):
    records=read_fasta(path)
    return [{"sequence_id":k,**sequence_statistics(v)} for k,v in records.items()]

def extract_interval(path, sequence_id, start, end):
    """Extract 1-based inclusive sequence coordinates from FASTA."""
    start,end=validate_interval(start,end); records=read_fasta(path)
    if sequence_id not in records: raise KeyError(f"Sequence {sequence_id!r} not found")
    seq=records[sequence_id]
    if end>len(seq): raise ValueError("Interval extends beyond sequence length")
    return seq[start-1:end]
