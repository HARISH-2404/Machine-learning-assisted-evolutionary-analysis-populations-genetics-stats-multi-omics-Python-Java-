"""DNA sequence validation, reverse-complement and composition."""
from collections import Counter

DNA = set("ACGTN")

def validate_dna(sequence: str, allow_n: bool = True) -> str:
    """Return uppercase DNA sequence or raise ValueError for invalid symbols."""
    seq = "".join(sequence.split()).upper().replace("U", "T")
    allowed = DNA if allow_n else set("ACGT")
    bad = sorted(set(seq) - allowed)
    if bad:
        raise ValueError(f"Invalid DNA symbols: {bad}")
    if not seq:
        raise ValueError("Sequence is empty")
    return seq

def reverse_complement(sequence: str) -> str:
    seq = validate_dna(sequence)
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def gc_content(sequence: str) -> float:
    seq = validate_dna(sequence)
    known = [b for b in seq if b != "N"]
    return 100.0 * sum(b in "GC" for b in known) / len(known) if known else float("nan")

def nucleotide_composition(sequence: str) -> dict:
    seq = validate_dna(sequence)
    counts = Counter(seq)
    return {base: counts.get(base, 0) for base in "ACGTN"}
