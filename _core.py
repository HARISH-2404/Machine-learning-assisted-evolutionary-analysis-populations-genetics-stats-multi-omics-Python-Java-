"""Shared helpers for FASTA, VCF and genomic interval operations."""
from pathlib import Path
import gzip
import numpy as np

def open_text(path):
    path=str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt", encoding="utf-8")

def read_fasta(path):
    """Read FASTA into {record_id: uppercase_sequence}; first header token is ID."""
    records={}; name=None; chunks=[]
    with open_text(path) as handle:
        for line in handle:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: records[name]="".join(chunks).upper()
                name=line[1:].split()[0]; chunks=[]
                if not name: raise ValueError("Empty FASTA identifier")
            else:
                if name is None: raise ValueError("Sequence encountered before FASTA header")
                chunks.append(line)
    if name is not None: records[name]="".join(chunks).upper()
    return records

def parse_info(info):
    out={}
    if info in ("", "."): return out
    for token in info.split(";"):
        if "=" in token:
            k,v=token.split("=",1); out[k]=v
        else: out[token]=True
    return out

def parse_vcf(path):
    """Yield VCF records as dictionaries. Supports plain or gzip VCF and simple biallelic SNPs."""
    with open_text(path) as handle:
        samples=[]
        for line in handle:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples=line.rstrip().split("\t")[9:]
                continue
            if not line.strip() or line.startswith("#"): continue
            f=line.rstrip().split("\t")
            if len(f)<8: raise ValueError(f"Malformed VCF row: {line[:100]}")
            yield {"chrom":f[0],"pos":int(f[1]),"id":f[2],"ref":f[3],"alt":f[4],
                   "qual":None if f[5]=="." else float(f[5]),"filter":f[6],
                   "info":parse_info(f[7]),"format":f[8] if len(f)>8 else "",
                   "samples":f[9:] if len(f)>9 else [],"sample_names":samples}

def dosage_from_gt(gt):
    """Return alternate-allele dosage for diploid GT; missing/non-diploid -> NaN."""
    import math
    if gt in (".","./.",".|.","","./.|."): return math.nan
    alleles=gt.replace("|","/").split("/")
    if len(alleles)!=2 or "." in alleles: return math.nan
    try: return float(sum(int(a) for a in alleles))
    except ValueError: return math.nan

def parse_genotypes(record):
    """Extract diploid dosage from GT field for each sample."""
    import math
    fmt=record["format"].split(":") if record["format"] else []
    if "GT" not in fmt: return [math.nan]*len(record["samples"])
    idx=fmt.index("GT"); result=[]
    for sample in record["samples"]:
        vals=sample.split(":")
        result.append(dosage_from_gt(vals[idx]) if idx<len(vals) else math.nan)
    return result

def validate_interval(start, end):
    if int(start)!=start or int(end)!=end or start<1 or end<start:
        raise ValueError("Intervals must be 1-based inclusive with 1 <= start <= end")
    return int(start),int(end)
