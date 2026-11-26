"""Finlay-Wilkinson style environmental response/stability summaries."""
def stability_parameters(records, genotype="genotype", environment="environment", trait="trait"):
    """Regress each genotype's performance on environmental means.

    Returns slope (bi), intercept and R2. Requires >=3 environments.
    """
    envvals={}
    for r in records: envvals.setdefault(r[environment],[]).append(float(r[trait]))
    envmean={e:sum(v)/len(v) for e,v in envvals.items()}
    if len(envmean)<3: raise ValueError("At least three environments are required")
    grand=sum(envmean.values())/len(envmean)
    x=[v-grand for v in envmean.values()]
    sx=sum(v*v for v in x)
    if sx==0: raise ValueError("No environmental variation")
    grouped={}
    for r in records: grouped.setdefault(r[genotype],{}).setdefault(r[environment],[]).append(float(r[trait]))
    out={}
    for g, byenv in grouped.items():
        pairs=[(x[i],sum(byenv[e])/len(byenv[e])) for i,e in enumerate(envmean) if e in byenv]
        if len(pairs)<3: continue
        mx=sum(a for a,b in pairs)/len(pairs); my=sum(b for a,b in pairs)/len(pairs)
        den=sum((a-mx)**2 for a,b in pairs)
        if den==0: continue
        slope=sum((a-mx)*(b-my) for a,b in pairs)/den
        intercept=my-slope*mx
        pred=[intercept+slope*a for a,b in pairs]
        ssr=sum((b-p)**2 for (a,b),p in zip(pairs,pred))
        sst=sum((b-my)**2 for a,b in pairs)
        out[g]={"intercept":intercept,"environmental_slope":slope,
                "r_squared":1-ssr/sst if sst else float("nan"),
                "mean_performance":my}
    return out
