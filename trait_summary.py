"""Descriptive statistics for numeric trait observations."""
import math
from collections import defaultdict

def trait_summary(records, trait, group_by=None):
    """Summarize trait values; records are dictionaries."""
    groups=defaultdict(list)
    for row in records:
        value=row.get(trait)
        if value is None: continue
        x=float(value)
        if not math.isfinite(x): continue
        key=row.get(group_by,"ALL") if group_by else "ALL"
        groups[key].append(x)
    out={}
    for key, vals in groups.items():
        n=len(vals); mean=sum(vals)/n
        var=sum((x-mean)**2 for x in vals)/(n-1) if n>1 else float("nan")
        out[str(key)]={"n":n,"mean":mean,"sd":math.sqrt(var) if n>1 else float("nan"),
                       "min":min(vals),"max":max(vals),"cv_percent":100*math.sqrt(var)/mean if n>1 and mean else float("nan")}
    return out
