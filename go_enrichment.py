"""Over-representation analysis for Gene Ontology term annotations.

Input mapping: gene ID -> iterable of GO IDs. Uses hypergeometric upper-tail
tests and Benjamini-Hochberg correction. Supply the tested/background universe.
"""
from collections import defaultdict
from ._enrichment import hypergeom_sf, bh_adjust

def go_enrichment(selected_genes, gene_to_go, background_genes=None,
                  min_term_size=2, max_term_size=None):
    selected=set(selected_genes)
    background=set(background_genes) if background_genes is not None else set(gene_to_go)
    selected &= background
    if not selected or not background:
        return []
    term_genes=defaultdict(set)
    for gene in background:
        for term in gene_to_go.get(gene, []):
            term_genes[str(term)].add(gene)
    M=len(background); N=len(selected); rows=[]
    for term, members in term_genes.items():
        n=len(members)
        if n<min_term_size or (max_term_size is not None and n>max_term_size): continue
        overlap=selected & members; k=len(overlap)
        if not k: continue
        p=hypergeom_sf(k-1,M,n,N)
        rows.append({"term":term,"overlap":k,"term_size":n,"selected_size":N,
                     "background_size":M,"genes":sorted(overlap),
                     "fold_enrichment":(k/N)/(n/M),"p_value":p})
    q=bh_adjust([r["p_value"] for r in rows])
    for row,adj in zip(rows,q): row["p_adjusted"]=adj
    return sorted(rows,key=lambda r:(r["p_adjusted"],r["p_value"],r["term"]))
