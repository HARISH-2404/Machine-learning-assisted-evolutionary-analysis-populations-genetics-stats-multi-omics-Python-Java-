"""Deterministic migration among demes."""
import numpy as np

def migration_step(frequencies,migration_matrix):
    """Rows of M are source demes, columns destination demes; each row sums to 1."""
    p=np.asarray(frequencies,float); m=np.asarray(migration_matrix,float)
    if p.ndim!=1 or m.shape!=(len(p),len(p)): raise ValueError('migration matrix must be square and match deme count')
    if np.any((p<0)|(p>1)) or np.any((m<0)|(m>1)) or not np.all(np.isfinite(m)): raise ValueError('values must be finite and in [0,1]')
    if not np.allclose(m.sum(axis=1),1): raise ValueError('each source row must sum to one')
    return p@m

def simulate_migration(initial_frequencies,migration_matrix,generations):
    if int(generations)!=generations or generations<0: raise ValueError('generations must be nonnegative integer')
    p=np.asarray(initial_frequencies,float); rows=[p.copy()]
    for _ in range(int(generations)): p=migration_step(p,migration_matrix); rows.append(p.copy())
    return np.vstack(rows)
