"""Join gene identifiers to functional annotations from a tabular file.

Expected TSV columns: gene_id, description, go_terms, kegg_pathways.
GO/KEGG term cells are semicolon-separated. This is a lightweight table
reader; it does not perform sequence homology or domain-based annotation.
"""
import csv
from pathlib import Path

def load_annotations(path, delimiter="\t"):
    required={"gene_id"}
    records={}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader=csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Annotation file must contain a gene_id column")
        for row in reader:
            gid=(row.get("gene_id") or "").strip()
            if not gid: continue
            records[gid]={
                "gene_id":gid,
                "description":(row.get("description") or "").strip(),
                "go_terms":[x.strip() for x in (row.get("go_terms") or "").split(";") if x.strip()],
                "kegg_pathways":[x.strip() for x in (row.get("kegg_pathways") or "").split(";") if x.strip()]
            }
    return records

def annotate_gene_list(gene_ids, annotations):
    result=[]
    for gid in gene_ids:
        row=annotations.get(gid, {"gene_id":gid,"description":"","go_terms":[],"kegg_pathways":[]})
        result.append(dict(row, annotation_status="annotated" if row["description"] or row["go_terms"] or row["kegg_pathways"] else "unannotated"))
    return result

def annotation_coverage(gene_ids, annotations):
    genes=set(gene_ids)
    hit=genes.intersection(annotations)
    return {"input_genes":len(genes),"annotated_genes":len(hit),
            "unannotated_genes":sorted(genes-hit),
            "coverage_percent":100*len(hit)/len(genes) if genes else 0.0}
