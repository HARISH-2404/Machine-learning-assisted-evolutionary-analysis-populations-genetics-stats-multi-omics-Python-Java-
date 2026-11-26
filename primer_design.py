"""Transparent heuristic primer metrics and simple candidate enumeration.

This is educational screening, not a replacement for Primer3 or in-silico PCR.
"""
from .sequence import validate_dna, reverse_complement, gc_content

def primer_metrics(primer: str) -> dict:
    seq = validate_dna(primer, allow_n=False)
    gc = gc_content(seq)
    # Wallace rule; approximate only for short oligos.
    tm = 2*(seq.count("A")+seq.count("T")) + 4*(seq.count("G")+seq.count("C"))
    return {"sequence":seq, "length":len(seq), "gc_percent":gc, "tm_wallace_c":tm,
            "reverse_complement":reverse_complement(seq),
            "homopolymer_max":max((len(run) for run in __import__("re").findall(r"(A+|C+|G+|T+)", seq)), default=0)}

def design_primers(template: str, min_len: int = 18, max_len: int = 22,
                   min_gc: float = 35, max_gc: float = 65,
                   min_tm: float = 50, max_tm: float = 65) -> dict:
    """Return first/last terminal candidates passing broad length/GC/Tm filters."""
    seq = validate_dna(template, allow_n=False)
    if min_len < 10 or max_len < min_len:
        raise ValueError("Use min_len >= 10 and max_len >= min_len")
    def candidates(region):
        out=[]
        for length in range(min_len, max_len+1):
            if len(region) < length: continue
            p=region[:length]
            m=primer_metrics(p)
            if min_gc <= m["gc_percent"] <= max_gc and min_tm <= m["tm_wallace_c"] <= max_tm and m["homopolymer_max"] <= 4:
                out.append(m)
        return out
    forward=candidates(seq)
    reverse=candidates(reverse_complement(seq))
    return {"forward_candidates":forward, "reverse_candidates":reverse,
            "note":"Candidates are not validated for dimers, hairpins, specificity, amplicon size or genome-wide uniqueness."}
