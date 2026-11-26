"""KEGG pathway over-representation analysis.

Provide a gene-to-pathway mapping prepared for the organism and KEGG release
being analysed. The module does not query KEGG or infer orthology.
"""
from .go_enrichment import go_enrichment

def kegg_enrichment(selected_genes, gene_to_pathway, background_genes=None,
                    min_pathway_size=2, max_pathway_size=None):
    rows=go_enrichment(selected_genes,gene_to_pathway,background_genes,
                       min_term_size=min_pathway_size,
                       max_term_size=max_pathway_size)
    for row in rows:
        row["pathway_id"]=row.pop("term")
    return rows
