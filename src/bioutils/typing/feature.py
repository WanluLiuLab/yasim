"""
feature.py -- General-Purposed GTF/GFF3/BED Record that Represents a Genomic Feature

This module includes GTF/GFF3/BED record datastructure and their one-line parsers.
"""

from __future__ import annotations

import uuid
from abc import abstractmethod
from typing import Union, Optional, Dict

from commonutils.stdlib_helper.logger_helper import get_logger
from commonutils.str_utils import to_dict

lh = get_logger(__name__)

GTFAttributeType = Dict[str, Union[str, int, float, bool, None]]
"""Type of GTF/GFF fields"""

GFF3_TOPLEVEL_NAME = "YASIM_GFF_TOPLEVEL"
"""The top-level virtual parent for GFF record that does not have a parent"""

VAILD_GTF_QUOTE_OPTONS = (
    "none",
    "blank",
    "string",
    "all"
)
"""
Valid GTF Quoting Options. They are:

* "none": Never quote, even if blanks found inside.
* "blank": Quote if blanks (\\r, \\n, \\t, space, etc.) found inside.
* "string": Quote if the field have type string. Will not quote for numeric types.
* "all": Quote all fields.
"""


class FeatureType(object):
    """
    Abstract type of a general GTF/GFF/BED Record.
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
        'attribute',
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
    Positive (``+``) or negative(``-``) or unknown (``.``)
    """

    frame: str
    """
    frame: One of ``0`` (first base of the feature is the first base of a codon),
                    ``1`` (the second base is the first base of a codon) or ``2``.
    """

    attribute: GTFAttributeType
    """Other attributes presented in Key-Value pair"""

    @abstractmethod
    def __eq__(self, other: FeatureType):
        pass

    @abstractmethod
    def __ne__(self, other: FeatureType):
        pass

    @abstractmethod
    def overlaps(self, other: FeatureType) -> bool:
        pass

    @abstractmethod
    def __gt__(self, other: FeatureType):
        pass

    @abstractmethod
    def __ge__(self, other: FeatureType):
        pass

    @abstractmethod
    def __lt__(self, other: FeatureType):
        pass

    @abstractmethod
    def __le__(self, other: FeatureType):
        pass

    @abstractmethod
    def format_string(self, **kwargs) -> str:
        pass


class Feature(FeatureType):
    """
    A general GTF/GFF/BED Record.
    """

    def __init__(self,
                 seqname: str,
                 source: str,
                 feature: str,
                 start: int,
                 end: int,
                 score: Union[int, float],
                 strand: str,
                 frame: str,
                 attribute: Optional[GTFAttributeType] = None
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
        if attribute is None:
            attribute = {}
        self.attribute = attribute

    def __eq__(self, other: Feature):
        return self.start == other.start and \
               self.end == other.end and \
               self.seqname == other.seqname \
               and self.strand == other.strand

    def __ne__(self, other: Feature):
        return not self == other

    def overlaps(self, other: Feature) -> bool:
        if self.seqname != other.seqname:
            return False
        return self.start < other.start < self.end or \
               self.start < other.end < self.end or \
               (
                       other.start < self.start and
                       self.end < other.end
               )

    def __gt__(self, other: Feature):
        return self.seqname > other.seqname or (
                self.seqname == other.seqname and self.start > other.start
        )

    def __ge__(self, other: Feature):
        return self > other or self == other

    def __lt__(self, other: Feature):
        return self.seqname < other.seqname or (
                self.seqname == other.seqname and self.start < other.start
        )

    def __le__(self, other: Feature):
        return self < other or self == other

    @classmethod
    def from_string(cls, in_str: str):
        """
        To generate ONE record from existing string.

        :param in_str: Input dictionary.
        """
        pass

    @abstractmethod
    def format_string(self, **kwargs) -> str:
        pass

    def __repr__(self):
        return self.format_string()


class Gff3Record(Feature):
    """
    A general GTF Record.
    """

    id: str
    """The ID field in Gff3, cannot be none. Will initialized to UUID if not present"""

    parent_id: str
    """
    FIXME: Multi parents
    
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    """

    def __init__(self,
                 seqname: str,
                 source: str,
                 feature: str,
                 start: int,
                 end: int,
                 score: float,
                 strand: str,
                 frame: str,
                 attribute: GTFAttributeType):
        super(Gff3Record, self).__init__(
            seqname=seqname,
            source=source,
            feature=feature,
            start=start,
            end=end,
            score=score,
            strand=strand,
            frame=frame,
            attribute=attribute
        )
        self.id = attribute.get("ID", uuid.uuid4())
        self.parent_id = attribute.get("Parent", GFF3_TOPLEVEL_NAME)

    @classmethod
    def from_string(cls, in_str: str):
        """
        To generate ONE record from existing string.

        :param in_str: Input dictionary.
        """
        global lh
        lh.debug(f'Adding {in_str}')
        line_split = in_str.split('\t')

        required_fields = line_split[0:-1]
        attributes = to_dict(line_split[-1], field_sep='=', record_sep=';', quotation_mark='\"\'', resolve_str=True)

        # Score should be an integer
        if required_fields[5] == ".":
            required_fields[5] = "0"
        return Gff3Record(seqname=required_fields[0],
                          source=required_fields[1],
                          feature=required_fields[2],
                          start=int(required_fields[3]),
                          end=int(required_fields[4]),
                          score=int(float(required_fields[5])),
                          strand=(required_fields[6]),
                          frame=(required_fields[7]),
                          attribute=attributes)

    def format_string(self):
        attribute_str = ""
        for k, v in self.attribute.items():
            v_str = repr(v).replace("'", '"')
            attribute_str = f"{attribute_str}{k}={v_str}; "
        return ("\t".join((
            self.seqname,
            self.source,
            self.feature,
            str(self.start),
            str(self.end),
            str(self.score),
            self.strand,
            self.frame,
            attribute_str
        )))


class GtfRecord(Feature):
    """
    A general GTF Record.
    
    >>> gtf_str = 'chr1\\thg38_rmsk\\texon\\t50331337\\t50332274\\t1587.000000\\t+\\t.\\tgene_id "HAL1"; transcript_id "HAL1"; '
    >>> gtf_from_line = GtfRecord.from_string(gtf_str)
    >>> gtf_from_line.seqname
    'chr1'
    >>> gtf_from_line.source
    'hg38_rmsk'
    >>> gtf_from_line.feature
    'exon'
    >>> gtf_from_line.start
    50331337
    >>> gtf_from_line.end
    50332274
    >>> gtf_from_line.score
    1587
    >>> gtf_from_line.strand
    '+'
    >>> gtf_from_line.attribute['gene_id']
    'HAL1'
    >>> str(gtf_from_line)
    'chr1\\thg38_rmsk\\texon\\t50331337\\t50332274\\t1587\\t+\\t.\\tgene_id "HAL1"; transcript_id "HAL1"; '
    """

    @classmethod
    def from_string(cls, in_str: str):
        global lh
        in_str = in_str.rstrip('\n\r')
        lh.trace(f'Adding {in_str}')
        line_split = in_str.split('\t')

        required_fields = line_split[0:-1]
        attributes = to_dict(line_split[-1], field_sep=' ', record_sep=';', quotation_mark='\"\'', resolve_str=True)

        # Score should be an integer
        if required_fields[5] == ".":
            required_fields[5] = "0"
        return GtfRecord(seqname=required_fields[0],
                         source=required_fields[1],
                         feature=required_fields[2],
                         start=int(required_fields[3]),
                         end=int(required_fields[4]),
                         score=int(float(required_fields[5])),
                         strand=(required_fields[6]),
                         frame=(required_fields[7]),
                         attribute=attributes)

    def format_string(
            self,
            quote: str = "all"
    ):
        if quote not in VAILD_GTF_QUOTE_OPTONS:
            raise ValueError(f"Invalid quoting option {quote}, should be one in {VAILD_GTF_QUOTE_OPTONS}.")
        attribute_full_str = ""
        for k, v in self.attribute.items():
            attr_str = repr(v).replace("'", '')
            if quote == "blank":
                if "\r" in attr_str or "\n" in attr_str or "\f" in attr_str or "\t" in attr_str or " " in attr_str:
                    attr_str = f"\"{attr_str}\""
            elif quote == "string":
                if not isinstance(v, int) and not isinstance(v, float):
                    attr_str = f"\"{attr_str}\""
            elif quote == "all":
                attr_str = f"\"{attr_str}\""
            attribute_full_str = f"{attribute_full_str}{k} " + attr_str + "; "
        return ("\t".join((
            self.seqname,
            self.source,
            self.feature,
            str(self.start),
            str(self.end),
            str(self.score),
            self.strand,
            self.frame,
            attribute_full_str
        )))
