"""Transcription and standard-code translation (5' to 3')."""
from .sequence import validate_dna

CODON_TABLE = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def transcribe(dna: str) -> str:
    return validate_dna(dna).replace("T", "U")

def translate(sequence: str, frame: int = 0, stop_at_stop: bool = False,
              strict: bool = False) -> str:
    """Translate DNA/RNA; frame is 0, 1 or 2. Ambiguous codons become X."""
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1, or 2")
    seq = "".join(sequence.split()).upper().replace("U", "T")
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError("Sequence must contain A/C/G/T/N only")
    protein = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "X")
        if strict and aa == "X":
            raise ValueError(f"Ambiguous codon at position {i+1}: {codon}")
        if aa == "*" and stop_at_stop:
            break
        protein.append(aa)
    return "".join(protein)
