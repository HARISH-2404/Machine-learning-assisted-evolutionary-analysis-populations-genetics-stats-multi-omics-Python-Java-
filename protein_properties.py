"""Basic protein sequence composition and approximate molecular mass."""
AA_MASS = {"A":89.09,"R":174.20,"N":132.12,"D":133.10,"C":121.15,
"E":147.13,"Q":146.15,"G":75.07,"H":155.16,"I":131.17,"L":131.17,
"K":146.19,"M":149.21,"F":165.19,"P":115.13,"S":105.09,"T":119.12,
"W":204.23,"Y":181.19,"V":117.15}
def protein_properties(sequence: str) -> dict:
    seq="".join(sequence.split()).upper()
    if not seq or set(seq)-set(AA_MASS):
        raise ValueError("Protein must contain the 20 standard amino-acid letters only")
    from collections import Counter
    counts=Counter(seq); n=len(seq)
    # Peptide mass: sum free amino-acid masses minus water for peptide bonds.
    mass=sum(AA_MASS[a]*c for a,c in counts.items()) - 18.015*(n-1)
    return {"length_aa":n,"molecular_mass_da_approx":mass,
            "amino_acid_counts":{a:counts.get(a,0) for a in sorted(AA_MASS)},
            "fractional_composition":{a:counts.get(a,0)/n for a in sorted(AA_MASS)},
            "charged_fraction":sum(counts[a] for a in "DEKRH")/n}
