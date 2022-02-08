from typing import Dict, Union

from commonutils import logger
from commonutils.str_utils import to_dict

lh = logger.get_logger(__name__)

__version__ = 0.1

cdef class GtfRecord(object):
    """
    A general GTF Record.
    """

    cdef public str seqname
    cdef public str _sequence
    """
    The cached sequence.
    """

    cdef public str source
    cdef public str feature
    cdef public int start
    cdef public int end
    cdef public float score
    cdef public str strand
    cdef public str frame
    cdef public dict attribute
    cdef public gtf_handler
    """
    Which Gtf* does this record belongs to?
    """

    def __init__(self,
                 str seqname,
                 str source,
                 str feature,
                 int start,
                 int end,
                 float score,
                 str strand,
                 str frame,
                 attribute: Dict[str, Union[str,int , float]]):
        """
        The filenames are named after Ensembl specifications.

        .. warning::
            Ensembl uses different way to represent 5'UTR.

        :param seqname: Chromosome or Contig name.
        :param source: The source of this record. e.g. ``hg38_rmsk`` or ``ensembl``.
        :param feature: Feature type name. e.g. ``exon`` or ``start_codon`` or ``5UTR``.
        :param start: Inclusive 1-based start position.
        :param end: Inclusive 2-based start position.
        :param score: Some kind of scoring.
        :param strand: Positive (``+``) or negative(``-``)
        :param frame: One of ``0`` (first base of the feature is the first base of a codon),
                    ``1`` (the second base is the first base of a codon) or ``2``.
        :param attribute: A semicolon-separated list of tag-value pairs,
                          providing additional information about each feature.
                          e.g. ``transcript_id`` or ``gene_name``.
        """
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = attribute
        self.gtf_handler = None

    property sequence:
        def __get__(self):
            """
            Get cached sequence, or get one from Fasta.
            """
            cdef int i = 0
            if self._sequence == None and self.gtf_handler and self.gtf_handler.attached_fasta:
                try:
                    self._sequence = self.gtf_handler.attached_fasta.sequence(self.seqname, self.start - 1, self.end)
                except ValueError:
                    pass
            return self._sequence

    @classmethod
    def from_string(cls, str in_str):
        """
        To generate ONE record from existing string.

        :param in_str: Input dictionary.
        """
        global lh
        lh.debug(f'Adding {in_str}')
        line_split = in_str.split('\t')

        required_fields = line_split[0:-1]
        attributes = to_dict(line_split[-1], field_sep=' ', record_sep=';', quotation_mark='\"\'', resolve_str=True)
        if required_fields[5] == ".":
            required_fields[5] = 0
        return GtfRecord(seqname=required_fields[0],
                         source=required_fields[1],
                         feature=required_fields[2],
                         start=int(required_fields[3]),
                         end=int(required_fields[4]),
                         score=int(float(required_fields[5])),
                         strand=(required_fields[6]),
                         frame=(required_fields[7]),
                         attribute=attributes)

    def __repr__(self):
        """
        Generate GTF string.

        :return: GTF string.
        """
        attribute_str = ""
        for k, v in self.attribute.items():
            attribute_str = f"{attribute_str}{k} " + repr(v).replace("'", '"') + "; "
        return ("\t".join((
            self.seqname,
            (self.source),
            self.feature,
            str(self.start),
            str(self.end),
            str(self.score),
            self.strand,
            self.frame,
            attribute_str
        )))

    __str__ = __repr__
