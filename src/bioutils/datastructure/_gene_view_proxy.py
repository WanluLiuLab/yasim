"""
_gene_view_proy -- Purposed GTF/GFF3/BED Proxy for Features in GeneView without Data Loss
"""

from __future__ import annotations

import sys
from abc import abstractmethod
from typing import List, Dict, Callable, Optional, Type

from bioutils.algorithm.sequence import complement, reverse_complement
from bioutils.typing.feature import GtfRecord, Feature, FeatureType, GTFAttributeType, Gff3Record


class _BaseFeature(FeatureType):
    __slots__ = (
        "_data"
    )

    _data: Feature

    @property
    def start(self) -> int:
        return self._data.start

    @start.setter
    def start(self, value: int):
        self._data.start = value

    @property
    def end(self) -> int:
        return self._data.end

    @end.setter
    def end(self, value: int):
        self._data.end = value

    @property
    def strand(self) -> str:
        return self._data.strand

    @strand.setter
    def strand(self, value: str):
        self._data.strand = value

    @property
    def source(self) -> str:
        return self._data.source

    @source.setter
    def source(self, value: str):
        self._data.source = value

    @property
    def seqname(self) -> str:
        return self._data.seqname

    @seqname.setter
    def seqname(self, value: str):
        self._data.seqname = value

    @property
    def frame(self) -> str:
        return self._data.frame

    @frame.setter
    def frame(self, value: str):
        self._data.frame = value

    @property
    def attribute(self) -> GTFAttributeType:
        return self._data.attribute

    @attribute.setter
    def attribute(self, value: GTFAttributeType):
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
    def _from_gtf(cls, feature: GtfRecord):
        new_instance = cls()
        new_instance._data = feature
        new_instance._setup()
        new_instance._setup_gtf()
        return new_instance

    @classmethod
    def _from_gff3(cls, feature: Gff3Record):
        new_instance = cls()
        new_instance._data = feature
        new_instance._setup()
        new_instance._setup_gff3()
        return new_instance

    @classmethod
    def from_feature(cls, feature: Feature):
        if isinstance(feature, GtfRecord):
            return cls._from_gtf(feature)
        elif isinstance(feature, Gff3Record):
            return cls._from_gff3(feature)
        else:
            raise NotImplementedError(f"Not implemented for {type(feature)}!")

    def __eq__(self, other: _BaseFeature):
        return self._data == other._data

    def __ne__(self, other: _BaseFeature):
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

    @property
    def transcript_id(self) -> str:
        return self._data.attribute["transcript_id"]

    @transcript_id.setter
    def transcript_id(self, value: str):
        self._data.attribute["transcript_id"] = value

    @property
    def gene_id(self) -> str:
        return self._data.attribute["gene_id"]

    @gene_id.setter
    def gene_id(self, value: str):
        self._data.attribute["gene_id"] = value

    @property
    def exon_number(self) -> int:
        return self._data.attribute["exon_number"]

    @exon_number.setter
    def exon_number(self, value: int):
        self._data.attribute["exon_number"] = value

    def _setup_gtf(self) -> None:
        if "transcript_id" not in self._data.attribute:
            self._data.attribute["transcript_id"] = "UNKNOWN"
        if "gene_id" not in self._data.attribute:
            self._data.attribute["gene_id"] = "UNKNOWN"
        if "exon_number" not in self._data.attribute:
            self._data.attribute["exon_number"] = 0

    def _setup_gff3(self) -> None:
        raise NotImplementedError

    def __repr__(self):
        return f"Exon {self.exon_number} of {self.transcript_id}"


class Transcript(_BaseFeature):
    __slots__ = (
        "exons",
        "_cdna_sequence"
    )
    exons: List[Exon]
    _cdna_sequence: Optional[str]

    @property
    def transcript_id(self) -> str:
        return self._data.attribute["transcript_id"]

    @transcript_id.setter
    def transcript_id(self, value: str):
        self._data.attribute["transcript_id"] = value

    @property
    def gene_id(self) -> str:
        return self._data.attribute["gene_id"]

    @gene_id.setter
    def gene_id(self, value: str):
        self._data.attribute["gene_id"] = value

    def _setup(self):
        self._cdna_sequence = None
        self.exons = []

    def _setup_gtf(self) -> None:
        if "transcript_id" not in self._data.attribute:
            self._data.attribute["transcript_id"] = "UNKNOWN"
        if "gene_id" not in self._data.attribute:
            self._data.attribute["gene_id"] = "UNKNOWN"

    def _setup_gff3(self) -> None:
        raise NotImplementedError

    def cdna_sequence(self, sequence_func: Callable[[str, int, int], str]) -> str:
        if self._cdna_sequence is not None:
            return self._cdna_sequence
        self._cdna_sequence = ""
        if self.strand == '-':
            for exon in sorted(self.exons):
                # print(self.seqname, exon.start - 1, exon.end, exon.exon_number)
                self._cdna_sequence += reverse_complement(sequence_func(self.seqname, exon.start - 1, exon.end))
            # print()
        else:
            for exon in self.exons:
                self._cdna_sequence += sequence_func(self.seqname, exon.start - 1, exon.end)
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

    def __repr__(self):
        return f"Transcript {self.transcript_id} of {self.gene_id}"


class Gene(_BaseFeature):
    __slots__ = (
        "transcripts"
    )
    transcripts: Dict[str, Transcript]

    @property
    def gene_id(self) -> str:
        return self._data.attribute["gene_id"]

    @gene_id.setter
    def gene_id(self, value: str):
        self._data.attribute["gene_id"] = value

    def _setup(self):
        self.transcripts = {}

    def _setup_gtf(self) -> None:
        if "gene_id" not in self._data.attribute:
            self._data.attribute["gene_id"] = "UNKNOWN"

    def _setup_gff3(self) -> None:
        raise NotImplementedError

    def __repr__(self):
        return f"Gene {self.gene_id}"
