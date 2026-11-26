"""Simple gene model validation and exon/intron summaries from coordinates."""
def summarize_gene(exons: list) -> dict:
    """Exons are 1-based inclusive (start,end) pairs on one chromosome."""
    if not exons:
        raise ValueError("At least one exon is required")
    ordered=sorted((int(s),int(e)) for s,e in exons)
    if any(s < 1 or e < s for s,e in ordered):
        raise ValueError("Invalid 1-based inclusive exon coordinates")
    for (_, prev_end),(start,_) in zip(ordered, ordered[1:]):
        if start <= prev_end:
            raise ValueError("Overlapping exons are not accepted")
    introns=[(ordered[i][1]+1, ordered[i+1][0]-1) for i in range(len(ordered)-1)]
    return {"exon_count":len(ordered),"exon_lengths":[e-s+1 for s,e in ordered],
            "intron_count":len(introns),"intron_intervals":introns,
            "spliced_length":sum(e-s+1 for s,e in ordered),
            "genomic_span":ordered[-1][1]-ordered[0][0]+1}
