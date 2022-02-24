from bioutils.io.fasta._fasta_view import FastaView

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
