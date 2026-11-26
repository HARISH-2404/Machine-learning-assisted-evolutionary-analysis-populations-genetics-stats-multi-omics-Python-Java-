"""Exact and IUPAC DNA motif scanning on a single strand."""
import re
IUPAC = {"A":"A","C":"C","G":"G","T":"T","R":"[AG]","Y":"[CT]",
         "S":"[GC]","W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]",
         "D":"[AGT]","H":"[ACT]","V":"[ACG]","N":"[ACGTN]"}
def scan_motif(sequence: str, motif: str, overlapping: bool = True) -> list:
    seq = "".join(sequence.split()).upper().replace("U","T")
    motif = motif.upper().replace("U","T")
    if not seq or not motif or set(seq)-set("ACGTN") or set(motif)-set(IUPAC):
        raise ValueError("Invalid DNA sequence or IUPAC motif")
    pattern = "".join(IUPAC[x] for x in motif)
    regex = re.compile(f"(?=({pattern}))" if overlapping else pattern)
    return [{"start":m.start()+1, "end":m.start()+len(m.group(1) if overlapping else m.group()),
             "match":m.group(1) if overlapping else m.group()} for m in regex.finditer(seq)]
