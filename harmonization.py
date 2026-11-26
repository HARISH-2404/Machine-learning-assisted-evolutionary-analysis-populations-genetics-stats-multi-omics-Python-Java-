"""Sample identifier harmonization and matching."""
import re
import pandas as pd

def harmonize_sample_ids(ids, lower=True, remove_spaces=True, strip_suffixes=()):
    """Normalize IDs deterministically; review collisions before using results."""
    out=[]
    for value in ids:
        s=str(value).strip()
        if remove_spaces: s=re.sub(r"\\s+", "", s)
        if lower: s=s.lower()
        for suffix in strip_suffixes:
            if s.endswith(suffix): s=s[:-len(suffix)]
        out.append(s)
    if len(set(out)) != len(out):
        raise ValueError("ID normalization created duplicate IDs; inspect collisions")
    return out

def intersect_samples(blocks, require_all=True):
    """Return sample IDs shared by named DataFrame blocks, preserving first-block order."""
    if not blocks: raise ValueError("At least one omics block is required")
    sets=[set(df.index.astype(str)) for df in blocks.values()]
    shared=set.intersection(*sets) if require_all else set.union(*sets)
    first=list(next(iter(blocks.values())).index.astype(str))
    return [x for x in first if x in shared]
