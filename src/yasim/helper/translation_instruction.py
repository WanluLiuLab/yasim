from __future__ import annotations

import enum
import json
import random
import uuid
from abc import abstractmethod

from labw_utils.bioutils.datastructure.gene_tree import GeneTreeInterface
from labw_utils.bioutils.datastructure.fasta_view import FastaViewType
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer, get_reader

from labw_utils.typing_importer import Tuple, List, Union, Mapping, Any, Dict, Optional
from yasim.helper import depth_io, depth
from yasim.helper.transposon import TransposonDatabase

DEFAULT_WEIGHT_TRANSCRIPT = 80
DEFAULT_WEIGHT_TE = 20
DEFAULT_WEIGHT_STOP = 100
DEFAULT_MINIMAL_SEQ_LEN = 250
DEFAULT_MINIMAL_TRANSPOSON_LEN = 20
DEFAULT_MINIMAL_TRANSCRIPT_LEN = 120


class SimpleSerializable:
    seq: str

    @classmethod
    @abstractmethod
    def from_dict(cls, d: Mapping[str, Any]):
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> Mapping[str, Any]:
        raise NotImplementedError

    def __eq__(self, others: SimpleSerializable) -> bool:
        if isinstance(others, self.__class__):
            raise TypeError
        return self.seq == others.seq

    def __hash__(self) -> int:
        return hash(self.seq)


class SimpleExon(SimpleSerializable):
    src_gene_id: str
    seq: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(**d)

    def to_dict(self):
        return {
            "type": "SimpleExon",
            "src_gene_id": self.src_gene_id,
            "seq": self.seq,
        }

    def __init__(
        self,
        *,
        src_gene_id: str,
        seq: str,
        **kwargs,
    ) -> None:
        _ = kwargs
        del kwargs
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
    depth: float
    _seq: Optional[str]

    def __init__(self, l: List[Union[SimpleExon, SimpleTE]], depth: float) -> None:
        self.l = l
        self.depth = depth
        self._seq = None

    def to_dict(self):
        return {
            "type": "SimpleTranscript",
            "l": {str(k): v.to_dict() for k, v in enumerate(self.l)},
            "d": self.depth,
        }

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(
            l=[
                SimpleExon.from_dict(v) if v["type"] == "SimpleExon" else SimpleTE.from_dict(v) for v in d["l"].values()
            ],
            depth=d["d"],
        )

    @property
    def seq(self) -> str:
        if self._seq is None:
            self._seq = "".join(it.seq for it in self.l)
        return self._seq


class TranslationInstructionState(enum.Enum):
    TRANSCRIPT = 0
    TRANSPOSON = 1
    STOP = 2


class TranslationInstructionStateAutomata:
    _weights: Tuple[float, float, float]
    _rdg: random.SystemRandom

    def __init__(self, weights: Tuple[float, float, float]):
        self._rdg = random.SystemRandom()
        self._weights = weights

    def draw(self) -> TranslationInstructionState:
        return self._rdg.choices(
            (
                TranslationInstructionState.TRANSCRIPT,
                TranslationInstructionState.TRANSPOSON,
                TranslationInstructionState.STOP,
            ),
            weights=self._weights,
            k=1,
        )[0]


class TranslationInstruction(SimpleSerializable):
    transcripts: Dict[str, SimpleTranscript]

    def __init__(self, transcripts: Dict[str, SimpleTranscript]) -> None:
        self.transcripts = transcripts

    def to_dict(self):
        return {k: v.to_dict() for k, v in self.transcripts.items()}

    def to_fasta(self, dst_fasta_path: str):
        with FastaWriter(dst_fasta_path) as faw:
            for k, v in self.transcripts.items():
                faw.write(FastaRecord(seq_id=k, sequence=v.seq))

    def to_json(self, dst_json_path: str):
        with get_writer(dst_json_path, is_binary=False) as w:
            json.dump(self.to_dict(), w, indent=4)

    def to_depth(self, dst_depth_path: str):
        depth_data = {k: v.depth for k, v in self.transcripts.items()}
        depth_io.write_depth(depth_data, dst_depth_path, feature_name="TRANSCRIPT_ID")

    @classmethod
    def from_dict(cls, d: Dict[str, Any]):
        return cls({k: SimpleTranscript.from_dict(v) for k, v in d.items()})

    @classmethod
    def from_json(cls, src_json_file_path: str):
        with get_reader(src_json_file_path, is_binary=False) as r:
            return cls.from_dict(json.load(r))

    @property
    def seq(self) -> str:
        raise ValueError

    @classmethod
    def generate(
        cls,
        *,
        n: int,
        tedb: TransposonDatabase,
        gt: GeneTreeInterface,
        fav: FastaViewType,
        mu: float = depth.DEFAULT_MU,
        weight_transcript: float = DEFAULT_WEIGHT_TRANSCRIPT,
        weight_transposon: float = DEFAULT_WEIGHT_TE,
        weight_stop: float = DEFAULT_WEIGHT_STOP,
        minimal_seq_len: int = DEFAULT_MINIMAL_SEQ_LEN,
        minimal_transposon_len: int = DEFAULT_MINIMAL_TRANSPOSON_LEN,
        minimal_transcript_len: int = DEFAULT_MINIMAL_TRANSCRIPT_LEN,
        high_cutoff_ratio: float = depth.DEFAULT_HIGH_CUTOFF_RATIO,
        low_cutoff: float = depth.DEFAULT_LOW_CUTOFF
    ):
        """

        :param n: Number of sequences.
        :return: Generated instance
        """
        rdg = random.SystemRandom()

        def autoclip(_seq: str, _min_len: int) -> str:
            while True:
                start = rdg.randint(0, len(seq) - 1)
                end = rdg.randint(start, len(seq))
                if end - start + 1 > _min_len:
                    return _seq[start:end]

        collapsed_transcripts = [gene.collapse_transcript(True) for gene in tqdm(gt.gene_values, "Collapsing...")]
        tisa = TranslationInstructionStateAutomata((weight_transcript, weight_transposon, weight_stop))
        final_simple_transcripts: Dict[str, SimpleTranscript] = {}
        while len(final_simple_transcripts) < n:
            new_transcript = SimpleTranscript(l=[], depth=0)
            while True:
                state = tisa.draw()
                if state == TranslationInstructionState.STOP:
                    break
                elif state == TranslationInstructionState.TRANSCRIPT:
                    transcript_to_use = rdg.choice(collapsed_transcripts)
                    seq = transcript_to_use.transcribe(fav.sequence, fav.legalize_region_best_effort)
                    if len(seq) < 2 * minimal_transcript_len:
                        continue
                    seq = autoclip(seq, minimal_transcript_len)
                    new_transcript.l.append(SimpleExon(src_gene_id=transcript_to_use.gene_id, seq=seq))
                elif state == TranslationInstructionState.TRANSPOSON:
                    src_te_name, seq = tedb.draw()
                    if len(seq) < 2 * minimal_transposon_len:
                        continue
                    seq = autoclip(seq, minimal_transposon_len)
                    new_transcript.l.append(SimpleTE(src_te_name=src_te_name, seq=seq))

            if len(new_transcript.seq) < minimal_seq_len:
                continue
            final_simple_transcripts["chimeric_transcript-" + str(uuid.uuid4())] = new_transcript
        for k, v in depth.simulate_gene_level_depth_gmm(
            gene_names=final_simple_transcripts.keys(),
            mu=mu,
            low_cutoff=low_cutoff,
            high_cutoff_ratio=high_cutoff_ratio
        ).items():
            final_simple_transcripts[k].depth = v
        return cls(final_simple_transcripts)
