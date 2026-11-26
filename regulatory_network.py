"""Build and summarize directed gene-regulatory candidate networks.

Edges are user-supplied evidence records (regulator, target, evidence, score).
A network edge is not proof of causal regulation.
"""
from collections import defaultdict, deque

def build_regulatory_network(edges, min_score=None):
    nodes=set(); clean=[]
    for edge in edges:
        source=str(edge["regulator"]); target=str(edge["target"])
        score=float(edge.get("score",1.0))
        if min_score is not None and score < min_score: continue
        nodes.update((source,target))
        clean.append({"regulator":source,"target":target,
                      "evidence":str(edge.get("evidence","unspecified")),
                      "score":score})
    return {"nodes":sorted(nodes),"edges":clean}

def network_summary(network):
    indeg=defaultdict(int); outdeg=defaultdict(int)
    for e in network["edges"]:
        outdeg[e["regulator"]]+=1; indeg[e["target"]]+=1
    nodes=network["nodes"]
    return {"node_count":len(nodes),"edge_count":len(network["edges"]),
            "in_degree":{n:indeg[n] for n in nodes},
            "out_degree":{n:outdeg[n] for n in nodes},
            "candidate_hubs":sorted(nodes,key=lambda n:outdeg[n],reverse=True)}

def reachable_targets(network, regulator):
    adjacency=defaultdict(list)
    for e in network["edges"]: adjacency[e["regulator"]].append(e["target"])
    seen=set(); queue=deque(adjacency[str(regulator)])
    while queue:
        node=queue.popleft()
        if node not in seen:
            seen.add(node); queue.extend(adjacency[node])
    return sorted(seen)
