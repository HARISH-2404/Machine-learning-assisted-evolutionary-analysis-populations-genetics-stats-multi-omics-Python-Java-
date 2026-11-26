# Molecular genetics package (Python)

Sequence-level utilities for reproducible, exploratory molecular-genetics work.

## Modules
- `sequence.py`: DNA validation, reverse complement, GC% and base composition.
- `translation.py`: transcription and standard genetic-code translation (frames 0–2).
- `codon_usage.py`: codon counts and RSCU.
- `motif_scanning.py`: overlapping/non-overlapping IUPAC motif scanning.
- `primer_design.py`: heuristic primer metrics and terminal candidate enumeration.
- `protein_properties.py`: amino-acid composition and approximate peptide mass.
- `gene_structure.py`: exon/intron coordinate summaries.
- `regulatory_motifs.py`: promoter-window extraction and motif scanning.

## Quick start
```python
from popgenml.molecular_genetics import reverse_complement, translate, gc_content
print(reverse_complement("ATGCCGTAA"))
print(translate("ATGGCTTAA"))
print(gc_content("ATGCGC"))
```

Run tests from this folder's parent:
```bash
python -m unittest discover -s tests -v
```

## Important scientific limits
These are transparent educational utilities, not validated clinical or breeding-grade pipelines. Translation uses the standard nuclear code; organellar or lineage-specific codes require an alternative table. Motif scans are sequence matches, not evidence of transcription-factor binding. Promoter extraction assumes a supplied genomic sequence and 1-based TSS; strand-aware extraction must be handled upstream for minus-strand genes. Primer candidates use the Wallace Tm approximation and do not assess hairpins, primer dimers, specificity, repeat regions, amplicon size, or homeolog cross-amplification. Use Primer3 plus in-silico PCR/BLAST for real primer design. Protein mass is approximate and ignores post-translational modifications. For cotton, distinguish A- and D-subgenome homeologs and use the correct transcript/CDS, strand, genetic code, and reference assembly.

## Functional annotation and pathway/network modules
- `functional_annotation.py`: load a gene-to-description/GO/KEGG TSV, annotate a list, and calculate coverage.
- `go_enrichment.py`: GO over-representation using a hypergeometric upper-tail test and BH FDR.
- `kegg_enrichment.py`: pathway over-representation using the same statistical method and supplied gene-to-pathway mapping.
- `regulatory_network.py`: construct directed regulator-target graphs, summarize degrees and list reachable downstream nodes.

The enrichment universe should be the genes that could have been selected in the experiment (for example, all genes tested in the differential-expression analysis), not automatically the whole genome. Use organism- and release-matched GO/KEGG mappings; GO parent propagation and ontology-aware tests (e.g. topGO) are not implemented. KEGG access/licensing and ID conversion are outside this offline module. Regulatory edges must carry provenance/evidence; correlation alone does not establish regulation.
