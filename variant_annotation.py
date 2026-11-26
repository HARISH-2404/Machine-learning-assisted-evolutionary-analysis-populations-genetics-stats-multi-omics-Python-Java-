"""Coordinate-based variant consequence annotation against GFF3 gene features."""
from ._core import parse_vcf

def read_gff_genes(path):
    """Read GFF3 gene records; coordinates remain 1-based inclusive."""
    genes=[]
    with open(path,encoding="utf-8") as h:
        for line in h:
            if not line.strip() or line.startswith("#"): continue
            f=line.rstrip().split("\t")
            if len(f)!=9 or f[2].lower()!="gene": continue
            attrs={}
            for item in f[8].split(";"):
                if "=" in item:
                    k,v=item.split("=",1); attrs[k]=v
            genes.append({"chrom":f[0],"start":int(f[3]),"end":int(f[4]),
                          "strand":f[6],"gene_id":attrs.get("ID",attrs.get("Name","unknown")),
                          "attributes":attrs})
    return genes

def annotate_variant(chrom, pos, genes, upstream=2000, downstream=500):
    """Assign nearest overlapping gene and coarse region class (not coding consequence)."""
    pos=int(pos); hits=[]
    for g in genes:
        if g["chrom"]!=chrom: continue
        if g["start"]<=pos<=g["end"]:
            region="genic"
            if pos==g["start"] or pos==g["end"]: region="gene_boundary"
            hits.append((0,g,region))
        else:
            dist=min(abs(pos-g["start"]),abs(pos-g["end"]))
            if dist<=upstream: hits.append((dist,g,"near_gene"))
    if not hits: return {"gene_id":None,"region":"intergenic","distance":None}
    hits.sort(key=lambda x:x[0]); d,g,r=hits[0]
    return {"gene_id":g["gene_id"],"region":r,"distance":int(d)}

def annotate_vcf(vcf_path, gff_path, upstream=2000):
    genes=read_gff_genes(gff_path); rows=[]
    for v in parse_vcf(vcf_path):
        a=annotate_variant(v["chrom"],v["pos"],genes,upstream=upstream)
        rows.append({"chrom":v["chrom"],"pos":v["pos"],"variant_id":v["id"],**a})
    return rows
