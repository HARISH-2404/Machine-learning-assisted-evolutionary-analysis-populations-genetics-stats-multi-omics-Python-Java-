"""Internal exact hypergeometric upper-tail test and BH adjustment."""
import math
def hypergeom_sf(k_minus_1, M, n, N):
    """P[X >= k], X~Hypergeom(M,n,N); inputs are population, successes, draws."""
    k=k_minus_1+1
    lo=max(k, N-(M-n)); hi=min(n,N)
    if M < 0 or n < 0 or N < 0 or n > M or N > M:
        raise ValueError("Invalid hypergeometric parameters")
    if lo>hi: return 0.0
    den=math.lgamma(M+1)-math.lgamma(N+1)-math.lgamma(M-N+1)
    total=0.0
    for x in range(lo,hi+1):
        lp=(math.lgamma(n+1)-math.lgamma(x+1)-math.lgamma(n-x+1)
            +math.lgamma(M-n+1)-math.lgamma(N-x+1)-math.lgamma(M-n-(N-x)+1)-den)
        total += math.exp(lp)
    return min(1.0,total)

def bh_adjust(pvalues):
    m=len(pvalues)
    order=sorted(range(m),key=lambda i:pvalues[i])
    q=[1.0]*m; running=1.0
    for rank_idx in range(m-1,-1,-1):
        i=order[rank_idx]; rank=rank_idx+1
        running=min(running,pvalues[i]*m/rank)
        q[i]=min(1.0,running)
    return q
