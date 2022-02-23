from __future__ import annotations

from abc import abstractmethod
from typing import List, Dict, Callable, Optional

from bioutils.datastructure.gff_gtf_record import GtfRecord
from bioutils.io import fasta


class SimpleData:
    __slots__ = (
        "start",
        "end",
        "strand",
        "seqname"
    )

    start: int
    end: int
    strand: str
    seqname: str

    def __eq__(self, other: SimpleData):
        return self.start == other.start and \
               self.end == other.end and \
               self.seqname == other.seqname \
               and self.strand == other.strand

    def __ne__(self, other: SimpleData):
        return not self == other

    @classmethod
    @abstractmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        pass

    @abstractmethod
    def to_gtf_record(self) -> GtfRecord:
        return GtfRecord(
            seqname=self.seqname,
            source="GeneView",
            feature="exon",
            start=self.start,
            end=self.end,
            score=0,
            strand=self.strand,
            frame='.',
            attribute={}
        )

    def overlaps(self, other: SimpleData) -> bool:
        if self.seqname != other.seqname:
            return False
        return self.start < other.start < self.end or \
               self.start < other.end < self.end or \
               (
                       other.start < self.start and
                       self.end < other.end
               )
    def __repr__(self):
        return "SimpleData"

    def __str__(self):
        return repr(self)


class Exon(SimpleData):
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

class Transcript(SimpleData):
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
            self._cdna_sequence = fasta.complement(self._cdna_sequence)
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

class Gene(SimpleData):
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
