"""Parse transcript/gene features from GFF3 and map transcript IDs to genes."""
def parse_gff_attributes(text):
    attrs={}
    for token in text.strip().split(";"):
        if "=" in token:
            k,v=token.split("=",1); attrs[k]=v
    return attrs

def read_transcript_annotation(gff_path):
    """Return rows for gene, mRNA/transcript, exon and CDS features."""
    rows=[]
    with open(gff_path,encoding="utf-8") as h:
        for line in h:
            if not line.strip() or line.startswith("#"): continue
            f=line.rstrip().split("\t")
            if len(f)!=9: continue
            attrs=parse_gff_attributes(f[8])
            rows.append({"seqid":f[0],"source":f[1],"feature":f[2],"start":int(f[3]),
                         "end":int(f[4]),"score":f[5],"strand":f[6],"phase":f[7],
                         "id":attrs.get("ID"),"parent":attrs.get("Parent"),
                         "gene_id":attrs.get("gene_id",attrs.get("gene",attrs.get("Parent"))),
                         "attributes":attrs})
    return rows

def transcript_gene_map(gff_path):
    rows=read_transcript_annotation(gff_path)
    return {r["id"]:r["parent"] for r in rows if r["feature"].lower() in ("mrna","transcript") and r["id"]}
