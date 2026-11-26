"""Codon counts and relative synonymous codon usage (RSCU)."""
from collections import Counter, defaultdict
from .translation import CODON_TABLE

def codon_usage(sequence: str) -> dict:
    seq = "".join(sequence.split()).upper().replace("U", "T")
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError("Sequence must contain A/C/G/T/N only")
    counts = Counter(seq[i:i+3] for i in range(0, len(seq)-2, 3)
                     if "N" not in seq[i:i+3])
    families = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        if aa != "*":
            families[aa].append(codon)
    rscu = {}
    for aa, codons in families.items():
        total = sum(counts[c] for c in codons)
        for c in codons:
            rscu[c] = (counts[c] / (total / len(codons))) if total else 0.0
    return {"codon_counts": dict(sorted(counts.items())),
            "rscu": dict(sorted(rscu.items())),
            "complete_codons": sum(counts.values()),
            "excluded_ambiguous_or_incomplete": (len(seq)//3)-sum(counts.values())}
