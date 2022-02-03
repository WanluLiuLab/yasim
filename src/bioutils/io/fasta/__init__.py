from bioutils.io.fasta.fasta_view import FastaView

_comp_trans = str.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')


def complement(seq: str) -> str:
    """
    Get reverse-complement of a sequence
    """
    return seq.translate(_comp_trans)


def reverse_complement(seq: str) -> str:
    """
    Get reverse-complement of a sequence
    """
    return complement(seq)[::-1]
