"""Promoter-window extraction and motif-hit enrichment counts."""
from .motif_scanning import scan_motif
def promoter_sequence(sequence: str, tss_1based: int, upstream: int=2000,
                      downstream: int=200) -> str:
    seq="".join(sequence.split()).upper()
    if set(seq)-set("ACGTN") or not 1 <= tss_1based <= len(seq):
        raise ValueError("Invalid sequence or TSS coordinate")
    start=max(1,tss_1based-upstream); end=min(len(seq),tss_1based+downstream-1)
    return seq[start-1:end]
def promoter_motif_hits(promoter: str, motifs: list) -> dict:
    return {motif:scan_motif(promoter,motif) for motif in motifs}
