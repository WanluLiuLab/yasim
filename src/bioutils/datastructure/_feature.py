from typing import Union, Callable, Optional

from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

__version__ = 0.1


class Feature(object):
    """
    A general GTF Record.
    """

    __slots__ = (
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        '_sequence',
        'sequence_func'
    )

    seqname: str
    """
    Chromosome or Contig name.
    """

    _sequence: Optional[str]
    """
    The cached sequence.
    """

    source: str
    """
    The source of this record. e.g. ``hg38_rmsk`` or ``ensembl``.
    """

    feature: str
    """
    Feature type name. e.g. ``exon`` or ``start_codon`` or ``5UTR``.
    """

    start: int
    """
    Inclusive 1-based start position.
    """

    end: int
    """
    Inclusive 1-based end position."""

    score: Union[int, float]
    """
    Some kind of scoring.
    """

    strand: str
    """
    Positive (``+``) or negative(``-``)
    """

    frame: str
    """
    frame: One of ``0`` (first base of the feature is the first base of a codon),
                    ``1`` (the second base is the first base of a codon) or ``2``.
    """

    sequence_func: Optional[Callable[[str, int, int], str]]
    """
    Function to call to get sequence, need to provide seqname, start and end, and return sequence in bytes.
    """

    def __init__(self,
                 seqname: str,
                 source: str,
                 feature: str,
                 start: int,
                 end: int,
                 score: Union[int, float],
                 strand: str,
                 frame: str
                 ):
        """
        The filenames are named after Ensembl specifications.

        .. warning::
            Ensembl uses different way to represent 5'UTR.
        """
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.sequence_func = None

    def sequence(self, sequence_func: Callable[[str, int, int], str]):
        """
        Get cached sequence, or get one from Fasta.
        """
        if self._sequence is None and self.sequence_func is not None:
            self._sequence = self.sequence_func(self.seqname, self.start - 1, self.end)
        return self._sequence
