"""
gene_view_proy -- GTF/GFF3/BED Record Proxy for Features in GeneView without Data Loss
"""

from __future__ import annotations

import bisect
import copy
import math
import uuid
from abc import abstractmethod
from typing import List, Dict, Callable, Optional, Iterable, Tuple, Type

from bioutils.algorithm.sequence import reverse_complement
from bioutils.typing.feature import GtfRecord, Feature, FeatureType, GTFAttributeType, Gff3Record
from commonutils.dynamic.hook_helper import hookable_decorator


def unknown_transcript_id() -> str:
    """Generate a new unknown transcript ID"""
    return 'unknown_transcript_id' + str(uuid.uuid4())


def unknown_gene_id() -> str:
    """Generate a new unknown gene ID"""
    return 'unknown_gene_id' + str(uuid.uuid4())


VALID_SORT_EXON_EXON_STRAND_POLICY = ("unstranded", "stranded", "none")
DEFAULT_SORT_EXON_EXON_STRAND_POLICY = "unstranded"


class GVPError(ValueError): pass
class ExonInATranscriptOnDifferentChromosomeError(GVPError): pass
class DuplicatedExonError(GVPError): pass
class ExonInATranscriptOnDifferentStrandError(GVPError): pass

@hookable_decorator
class BaseFeatureProxy(FeatureType):
    """
    Base class of Feature Proxy.
    """

    _data: Feature

    _was_modified: bool

    def duplicate_cast(self, class_type: Type[BaseFeatureProxy]):
        return class_type.from_feature(copy.deepcopy(self._data))

    def copy_data(self):
        """
        Make a copy of data
        """
        self._data = copy.deepcopy(self._data)

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
    def feature(self) -> str:
        return self._data.feature

    @feature.setter
    def feature(self, value: str):
        self._data.feature = value

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

    def __eq__(self, other: BaseFeatureProxy):
        return self._data == other._data

    def __ne__(self, other: BaseFeatureProxy):
        return self._data != other._data

    def overlaps(self, other: BaseFeatureProxy) -> bool:
        return self._data.overlaps(other._data)

    def __gt__(self, other: BaseFeatureProxy):
        return self._data > other._data

    def __ge__(self, other: BaseFeatureProxy):
        return self._data >= other._data

    def __lt__(self, other: BaseFeatureProxy):
        return self._data < other._data

    def __le__(self, other: BaseFeatureProxy):
        return self._data <= other._data

    def __repr__(self):
        return "BaseFeatureProxy"

    def __str__(self):
        return repr(self)

    def format_string(self, **kwargs) -> str:
        """Disabled"""
        pass


class Exon(BaseFeatureProxy):

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
            self._data.attribute["transcript_id"] = unknown_transcript_id()
        if "gene_id" not in self._data.attribute:
            self._data.attribute["gene_id"] = unknown_gene_id()
        if "exon_number" not in self._data.attribute:
            self._data.attribute["exon_number"] = 0

    def _setup_gff3(self) -> None:
        raise NotImplementedError

    def __repr__(self):
        return f"Exon {self.exon_number} of {self.transcript_id}"


class Transcript(BaseFeatureProxy):
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
            self._data.attribute["transcript_id"] = unknown_transcript_id()
        if "gene_id" not in self._data.attribute:
            self._data.attribute["gene_id"] = unknown_gene_id()

    def _setup_gff3(self) -> None:
        raise NotImplementedError

    @property
    def naive_length(self) -> int:
        """
        The length on GTF
        """
        return self.end - self.start + 1

    @property
    def span_length(self) -> int:
        """
        The spanning length of all exons
        """
        exon_s_min = math.inf
        exon_e_max = - math.inf
        for exon in self.exons:
            exon_s_min = min(exon_s_min, exon.start)
            exon_e_max = max(exon_e_max, exon.end)
        return exon_e_max - exon_s_min + 1

    @property
    def transcribed_length(self) -> int:
        """
        Length after transcribed to cDNA
        """
        reti = 0
        for exon in self.exons:
            reti += exon.end - exon.start + 1
        return reti

    def cdna_sequence(self, sequence_func: Callable[[str, int, int], str]) -> str:
        if self._cdna_sequence is not None:
            return self._cdna_sequence
        self._cdna_sequence = ""
        if self.strand == '-':
            for exon in sorted(self.exons)[::-1]:
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
        return True

    def sort_exons(
            self,
            exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY,
            remove_exon_duplicates: bool = True
    ):
        if len(self.exons) == 0:
            return
        self.exons = sorted(self.exons)
        if remove_exon_duplicates and len(self.exons) >= 2:
            exon_id = 1
            while exon_id < len(self.exons):
                if self.exons[exon_id] == self.exons[exon_id -1]:
                    self.exons.pop(exon_id)
                else:
                    exon_id += 1
        if exon_number_policy == "stranded":
            if self.strand == '+':
                for i in range(len(self.exons)):
                    self.exons[i].exon_number = i + 1
            elif self.strand == '-':
                for i in range(len(self.exons)):
                    self.exons[len(self.exons) - i - 1].exon_number = i + 1
        elif exon_number_policy == "unstranded":
            for i in range(len(self.exons)):
                self.exons[i].exon_number = i + 1

    def fast_add_exon(
            self,
            exon: Exon
    ):
        self.add_exon(exon, False,False, False)

    def add_exon(
            self,
            exon:Exon,
            check_duplicate:bool=True,
            check_same_chrome:bool=True,
            check_same_strand:bool=True
    ):
        new_pos = bisect.bisect_left(self.exons, exon)
        if check_duplicate and new_pos < len(self.exons) and self.exons[new_pos] == self.exons[new_pos + 1]:
            raise DuplicatedExonError
        if check_same_chrome and exon.seqname != self.seqname:
            raise ExonInATranscriptOnDifferentChromosomeError
        if check_same_strand and exon.strand != self.strand and exon.strand != "." and self.strand != ".":
            raise ExonInATranscriptOnDifferentStrandError
        self.exons.insert(new_pos, exon)

    @property
    def exon_boundaries(self) -> Iterable[Tuple[int, int]]:
        for exon in self.exons:
            yield exon.start, exon.end

    @property
    def splice_sites(self) -> Iterable[Tuple[int, int]]:
        le = len(self.exons)
        for i in range(le - 1):
            yield self.exons[i].end, self.exons[i + 1].start

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return f"Transcript {self.transcript_id} of {self.gene_id}"


class Gene(BaseFeatureProxy):
    transcripts: Dict[str, Transcript]

    # TODO: Unit test needed

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
            self._data.attribute["gene_id"] = unknown_gene_id()

    def _setup_gff3(self) -> None:
        raise NotImplementedError

    def __repr__(self):
        return f"Gene {self.gene_id}"
