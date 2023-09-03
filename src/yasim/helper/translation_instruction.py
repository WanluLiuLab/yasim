from abc import abstractmethod
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord

from labw_utils.typing_importer import Tuple, List, Union, Mapping, Any, Dict, Optional

SimpleExon: Tuple[str, bool, int, int]
"""Contig ID, strand, start, end"""


class SimpleSerializable:
    @classmethod
    @abstractmethod
    def from_dict(cls, d: Mapping[str, Any]):
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> Mapping[str, Any]:
        raise NotImplementedError

    @property
    @abstractmethod
    def seq(self) -> str:
        raise NotImplementedError
    
    def __eq__(self, others: object) -> bool:
        if not hasattr(others, "seq"):
            raise TypeError
        return self.seq == others.seq
    
    def __hash__(self) -> int:
        return hash(self.seq)


class SimpleExon(SimpleSerializable):
    contig_name: str
    strand: bool
    start: int
    end: int
    src_gene_id: str
    seq: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(**d)

    def to_dict(self):
        return {
            "type": "SimpleExon",
            "contig_name": self.contig_name,
            "strand": self.strand,
            "start": self.start,
            "end": self.end,
            "src_gene_id": self.src_gene_id,
            "seq": self.seq,
        }

    def __init__(
        self,
        *,
        contig_name: str,
        strand: bool,
        start: int,
        end: int,
        src_gene_id: str,
        seq: str,
        **kwargs,
    ) -> None:
        _ = kwargs
        del kwargs
        self.contig_name = contig_name
        self.strand = strand
        self.start = start
        self.end = end
        self.src_gene_id = src_gene_id
        self.seq = seq


class SimpleTE(SimpleSerializable):
    src_te_name: str
    seq: str

    def __init__(
        self,
        src_te_name: str,
        seq: str,
        **kwargs,
    ):
        _ = kwargs
        del kwargs
        self.src_te_name = src_te_name
        self.seq = seq

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(**d)

    def to_dict(self):
        return {
            "type": "SimpleTE",
            "src_te_name": self.src_te_name,
            "seq": self.seq,
        }


class SimpleTranscript(SimpleSerializable):
    l: List[Union[SimpleExon, SimpleTE]]
    _seq: Optional[str]

    def __init__(self, l: List[Union[SimpleExon, SimpleTE]]) -> None:
        self.l = l
        self._seq = None

    def to_dict(self):
        return {str(k): v for k, v in enumerate(self.l)}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(
            [
                SimpleExon.from_dict(v)
                if v["type"] == "SimpleExon"
                else SimpleTE.from_dict(v)
                for v in d.values()
            ]
        )

    @property
    def seq(self) -> str:
        if self._seq is None:
            self._seq =  "".join(it.seq for it in self.l)
        return self._seq
        


class TranslationInstruction(SimpleSerializable):
    transcripts: Dict[str, SimpleTranscript]

    def __init__(self, transcripts: Dict[str, SimpleTranscript]) -> None:
        self.transcripts = transcripts

    def to_dict(self):
        return self.transcripts

    def to_fasta(self, dst_fasta_path: str):
        with FastaWriter(dst_fasta_path) as faw:
            for k, v in self.transcripts.items():
                faw.write(FastaRecord(seq_id=k, sequence=v.seq))

    @classmethod
    def from_dict(cls, d: Dict[str, Any]):
        return cls({k: SimpleTranscript.from_dict(v) for k, v in d.items()})

    @property
    def seq(self) -> str:
        raise ValueError
