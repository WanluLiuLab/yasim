from __future__ import annotations

from typing import List, Dict

from bioutils.io.gtf import GtfRecord

class SimpleData:
    start: int
    end: int
    strand: str = "."
    seqname: str

    def __eq__(self, other:SimpleData):
        return self.start == other.start and \
               self.end == other.end and \
               self.seqname == other.seqname \
               and self.strand == other.strand



class SimpleExon(SimpleData):


    @classmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        new_instance = cls()
        new_instance.start = gtf_record.start
        new_instance.end = gtf_record.end
        new_instance.seqname = gtf_record.seqname
        new_instance.strand = gtf_record.strand
        return new_instance



    def __ne__(self, other):
        return not self == other


class Transcript(SimpleData):
    exons: List[SimpleExon]
    name: str
    gene_name: str

    @classmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        new_instance = cls()
        new_instance.start = gtf_record.start
        new_instance.end = gtf_record.end
        try:
            transcript_id = gtf_record.attribute['transcript_id']
            new_instance.name = transcript_id
        except KeyError:
            new_instance.name = None
        try:
            gene_id = gtf_record.attribute['gene_id']
            new_instance.gene_name = gene_id
        except KeyError:
            new_instance.gene_name = None
        new_instance.exons = []
        new_instance.strand = gtf_record.strand
        new_instance.seqname = gtf_record.seqname
        return new_instance

    def transcribe_cdna(self, fasta_handler: FastaView) -> str:
        rets = ""
        for exon in self.exons:
            try:
                rets += fasta_handler.sequence(self.seqname, exon.start - 1, exon.end)
            except ValueError:
                pass
        if self.strand == '-':
            rets = fasta.complement(rets)
        return rets

    def __eq__(self, other):
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


class Gene:
    name: str
    start: int
    end: int
    transcripts: Dict[str, Transcript]
    strand: str
    seqname: str

    def __init__(self):
        pass

    @classmethod
    def from_gtf_record(cls, gtf_record: GtfRecord):
        new_instance = cls()
        new_instance.start = gtf_record.start
        new_instance.end = gtf_record.end
        try:
            gene_id = gtf_record.attribute['gene_id']
            new_instance.name = gene_id
        except KeyError:
            new_instance.name = None
        new_instance.strand = gtf_record.strand
        new_instance.transcripts = {}
        new_instance.seqname = gtf_record.seqname
        return new_instance
