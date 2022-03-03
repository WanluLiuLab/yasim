"""
_gene_view_proy -- Purposed GTF/GFF3/BED Proxy for Features in GeneView without Data Loss
"""



from __future__ import annotations

from abc import abstractmethod
from typing import List, Dict, Callable, Optional, Type

from bioutils.algorithm.sequence import complement
from bioutils.typing.feature import GtfRecord, Feature, FeatureType, GTFAttributeType, Gff3Record


class _BaseFeature(FeatureType):
    __slots__ = (
        "start",
        "end",
        "strand",
        "source",
        "seqname",
        "feature",
        "strand",
        "attribute",
        "_data"
    )

    _data:Feature

    @property
    def start(self) -> int:
        return self._data.start

    @start.setter
    def start(self, value:int):
        self._data.start = value

    @property
    def end(self) -> int:
        return self._data.end

    @end.setter
    def end(self, value:int):
        self._data.end = value

    @property
    def strand(self) -> str:
        return self._data.strand

    @strand.setter
    def strand(self, value:str):
        self._data.strand = value

    @property
    def source(self) -> str:
        return self._data.source

    @source.setter
    def source(self, value:str):
        self._data.source = value

    @property
    def seqname(self) -> str:
        return self._data.seqname

    @seqname.setter
    def seqname(self, value:str):
        self._data.seqname = value

    @property
    def frame(self) -> str:
        return self._data.frame

    @frame.setter
    def frame(self, value:str):
        self._data.frame = value

    @property
    def attribute(self) -> GTFAttributeType:
        return self._data.attribute

    @attribute.setter
    def attribute(self, value:GTFAttributeType):
        self._data.attribute = value

    def get_data(self) -> Feature:
        """Read-only _data"""
        return self._data

    @abstractmethod
    def _setup_gtf(self) -> None:
        """
        This method prepares underlying method to set the record up for GTF
        """
        pass

    @abstractmethod
    def _setup_gff3(self) -> None:
        """
        This method prepares underlying method to set the record up for GFF3
        """
        pass

    def _setup(self):
        """
        This method prepares the object for arbitrary :py:class:`Feature` types from `data`.
        """
        pass

    @classmethod
    def _from_gtf(cls, feature:GtfRecord):
        new_instance = cls()
        new_instance._data = feature
        new_instance._setup()
        new_instance._setup_gtf()
        return new_instance

    @classmethod
    def _from_gff3(cls, feature:Gff3Record):
        new_instance = cls()
        new_instance._data = feature
        new_instance._setup()
        new_instance._setup_gff3()
        return new_instance

    @classmethod
    def from_feature(cls, feature:Feature):
        if isinstance(feature, GtfRecord):
            return cls._from_gtf(feature)
        elif isinstance(feature, Gff3Record):
            return cls._from_gff3(feature)
        else:
            pass
        pass


    def __eq__(self, other: _BaseFeature):
        return self._data == other._data

    def __ne__(self, other:_BaseFeature):
        return self._data != other._data

    def overlaps(self, other: _BaseFeature) -> bool:
        return self._data.overlaps(other._data)

    def __gt__(self, other: _BaseFeature):
        return self._data > other._data

    def __ge__(self, other: _BaseFeature):
        return self._data >= other._data

    def __lt__(self, other: _BaseFeature):
        return self._data < other._data

    def __le__(self, other: _BaseFeature):
        return self._data <= other._data

    def __repr__(self):
        return "_BaseFeature"

    def __str__(self):
        return repr(self)


class Exon(_BaseFeature):
    transcript_id: str
    gene_id: str
    exon_number: int

    __slots__ = (
        "transcript_id",
        "gene_id",
        "exon_number"
    )

    @classmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        new_instance = cls()
        new_instance.start = gtf_record.start
        new_instance.end = gtf_record.end
        new_instance.transcript_id = gtf_record.attribute['transcript_id']
        new_instance.gene_id = gtf_record.attribute['gene_id']
        try:
            new_instance.exon_number = gtf_record.attribute['exon_number']
        except KeyError:
            new_instance.exon_number = 0
        new_instance.strand = gtf_record.strand
        new_instance.seqname = gtf_record.seqname
        return new_instance

    def to_gtf_record(self):
        return GtfRecord(
            seqname=self.seqname,
            source="GeneView",
            feature="exon",
            start=self.start,
            end=self.end,
            score=0,
            strand=self.strand,
            frame='.',
            attribute={
                "gene_id": self.gene_id,
                "transcript_id": self.transcript_id,
                "exon_number": self.exon_number
            }
        )

    def __repr__(self):
        return f"Exon {self.exon_number} of {self.transcript_id}"


class Transcript(_BaseFeature):
    __slots__ = (
        "exons",
        "transcript_id",
        "gene_id",
        "_cdna_sequence"
    )
    exons: List[Exon]
    transcript_id: str
    gene_id: str
    _cdna_sequence: Optional[str]

    @classmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        new_instance = cls()
        new_instance.start = gtf_record.start
        new_instance.end = gtf_record.end
        new_instance.transcript_id = gtf_record.attribute['transcript_id']
        new_instance.gene_id = gtf_record.attribute['gene_id']
        new_instance.exons = []
        new_instance.strand = gtf_record.strand
        new_instance.seqname = gtf_record.seqname
        new_instance._cdna_sequence = None
        return new_instance

    def cdna_sequence(self, sequence_func: Callable[[str, int, int], str]) -> str:
        if self._cdna_sequence is not None:
            return self._cdna_sequence
        self._cdna_sequence = ""
        for exon in self.exons:
            try:
                self._cdna_sequence += sequence_func(self.seqname, exon.start - 1, exon.end)
            except ValueError:
                pass
        if self.strand == '-':
            self._cdna_sequence = complement(self._cdna_sequence)
        return self._cdna_sequence

    def __eq__(self, other: Transcript):
        for exon_s, exon_o in zip(self.exons, other.exons):
            if not exon_s == exon_o:
                return False
        return True  # Do not require span length equality
        # return self.start == other.start and \
        #        self.end == other.end and \
        #        self.seqname == other.seqname \
        #        and self.strand == other.strand

    def __ne__(self, other):
        return not self == other

    def to_gtf_record(self):
        return GtfRecord(
            seqname=self.seqname,
            source="GeneView",
            feature="transcript",
            start=self.start,
            end=self.end,
            score=0,
            strand=self.strand,
            frame='.',
            attribute={
                "gene_id": self.gene_id,
                "transcript_id": self.transcript_id
            }
        )

    def __repr__(self):
        return f"Transcript {self.transcript_id} of {self.gene_id}"


class Gene(_BaseFeature):
    __slots__ = (
        "transcripts",
        "gene_id"
    )
    gene_id: str
    transcripts: Dict[str, Transcript]

    @classmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        new_instance = cls()
        new_instance.start = gtf_record.start
        new_instance.end = gtf_record.end
        new_instance.gene_id = gtf_record.attribute['gene_id']
        new_instance.transcripts = {}
        new_instance.strand = gtf_record.strand
        new_instance.seqname = gtf_record.seqname
        return new_instance

    def to_gtf_record(self):
        return GtfRecord(
            seqname=self.seqname,
            source="GeneView",
            feature="gene",
            start=self.start,
            end=self.end,
            score=0,
            strand=self.strand,
            frame='.',
            attribute={
                "gene_id": self.gene_id
            }
        )

    def __repr__(self):
        return f"Gene {self.gene_id}"
