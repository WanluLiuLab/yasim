_comp_trans = str.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')


def complement(seq: str) -> str:
    """
    Get reverse-complement of a sequence

    >>> complement("CTGACTGA")
    'GACTGACT'
    """
    return seq.translate(_comp_trans)


def reverse_complement(seq: str) -> str:
    """
    Get reverse-complement of a sequence

    >>> reverse_complement("CTGACTGA")
    'TCAGTCAG'
    """
    return complement(seq)[::-1]

def get_gc_percent(seq:str) -> float:
    if len(seq) == 0:
        return 0
    gc = 0
    for base in seq:
        if base in ("C", "G", "c", "g"):
            gc += 1
    return gc/len(seq)
