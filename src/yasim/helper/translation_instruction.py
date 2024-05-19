from __future__ import annotations

import enum
import json
import random
from abc import abstractmethod

from labw_utils.bioutils.datastructure.fasta_view import (
    FastaViewType,
    normalize_nt_sequence,
)
from labw_utils.bioutils.datastructure.gene_tree import GeneTreeInterface
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.datastructure.quantification_optimized_feature_index import (
    QuantificationOptimizedFeatureIndex,
)
from labw_utils.bioutils.datastructure.transposon import TransposonDatabase
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.bioutils.record.feature import strand_repr
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer, get_reader
from labw_utils.typing_importer import List, Union, Mapping, Any, Dict, Optional
from yasim.helper import depth_io, depth

DEFAULT_MINIMAL_SEQ_LEN = 250
DEFAULT_MINIMAL_TRANSPOSON_LEN = 20
DEFAULT_MINIMAL_TRANSCRIPT_LEN = 120


class FusionTypes(enum.Enum):
    FivePrimeFusion = 0
    ThreePrimeFusion = 1
    GeneGeneFusion = 2
    GeneOnly = 3
    TransposonOnly = 4


DEFAULT_WEIGHTS = {
    FusionTypes.FivePrimeFusion: 100,
    FusionTypes.ThreePrimeFusion: 100,
    FusionTypes.GeneGeneFusion: 0,
    FusionTypes.GeneOnly: 0,
    FusionTypes.TransposonOnly: 0,
}


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


class SimpleCoordinate(SimpleSerializable):
    contig: str
    strand: str
    start: int
    end: int

    def __init__(self, contig: str, strand: str, start: int, end: int):
        self.contig = contig
        self.strand = strand
        self.start = start
        self.end = end

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(**d)

    def to_dict(self):
        return {
            "contig": self.contig,
            "strand": self.strand,
            "start": self.start,
            "end": self.end,
        }


class SimpleExon(SimpleSerializable):
    src_gene_id: str
    seq: str
    coordinate: SimpleCoordinate

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(**d)

    def to_dict(self):
        return {
            "type": "SimpleExon",
            "src_gene_id": self.src_gene_id,
            "seq": self.seq,
            "coordinate": self.coordinate.to_dict(),
        }

    def __init__(
        self,
        *,
        src_gene_id: str,
        seq: str,
        coordinate: Union[SimpleCoordinate, Dict],
        **kwargs,
    ) -> None:
        _ = kwargs
        del kwargs
        self.src_gene_id = src_gene_id
        self.seq = seq
        self.coordinate = (
            coordinate
            if isinstance(coordinate, SimpleCoordinate)
            else SimpleCoordinate.from_dict(coordinate)
        )


class SimpleTE(SimpleSerializable):
    src_te_name: str
    seq: str
    coordinate: SimpleCoordinate

    def __init__(
        self,
        src_te_name: str,
        seq: str,
        coordinate: Union[SimpleCoordinate, Dict],
        **kwargs,
    ):
        _ = kwargs
        del kwargs
        self.src_te_name = src_te_name
        self.seq = seq
        self.coordinate = (
            coordinate
            if isinstance(coordinate, SimpleCoordinate)
            else SimpleCoordinate.from_dict(coordinate)
        )

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(**d)

    def to_dict(self):
        return {
            "type": "SimpleTE",
            "src_te_name": self.src_te_name,
            "seq": self.seq,
            "coordinate": self.coordinate.to_dict(),
        }


class SimpleTranscript(SimpleSerializable):
    l: List[Union[SimpleExon, SimpleTE]]
    d: float
    _seq: Optional[str]

    def __init__(self, l: List[Union[SimpleExon, SimpleTE]], d: float) -> None:
        self.l = l
        self.d = d
        self._seq = None

    def to_dict(self):
        return {
            "type": "SimpleTranscript",
            "l": {str(k): v.to_dict() for k, v in enumerate(self.l)},
            "d": self.d,
        }

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]):
        return cls(
            l=[
                (
                    SimpleExon.from_dict(v)
                    if v["type"] == "SimpleExon"
                    else SimpleTE.from_dict(v)
                )
                for v in d["l"].values()
            ],
            d=d["d"],
        )

    @property
    def seq(self) -> str:
        if self._seq is None:
            self._seq = "".join(it.seq for it in self.l)
        return self._seq


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
        depth_data = {k: v.d for k, v in self.transcripts.items()}
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
        tedb: Optional[TransposonDatabase],
        transposon_fi: QuantificationOptimizedFeatureIndex,
        transposon_gt: GeneTreeInterface,
        gt: GeneTreeInterface,
        fav: FastaViewType,
        mu: float = depth.DEFAULT_MU,
        disable_gmm: bool = False,
        minimal_seq_len: int = DEFAULT_MINIMAL_SEQ_LEN,
        minimal_transposon_len: int = DEFAULT_MINIMAL_TRANSPOSON_LEN,
        minimal_transcript_len: int = DEFAULT_MINIMAL_TRANSCRIPT_LEN,
        high_cutoff_ratio: float = depth.DEFAULT_HIGH_CUTOFF_RATIO,
        low_cutoff: float = depth.DEFAULT_LOW_CUTOFF,
    ):
        rdg = random.SystemRandom()

        final_simple_transcripts: Dict[str, SimpleTranscript] = {}
        pbar = tqdm(desc="Generating sequences...", total=n)
        choices = list(DEFAULT_WEIGHTS.keys())
        weights = list(DEFAULT_WEIGHTS.values())

        def add_transcript(_transcript_to_use: Transcript) -> Optional[SimpleExon]:
            seq_to_add = _transcript_to_use.transcribe(
                fav.sequence,
                fav.legalize_region_best_effort,
            )
            if len(seq_to_add) < minimal_transcript_len:
                return None
            else:
                return SimpleExon(
                    src_gene_id=_transcript_to_use.gene_id,
                    seq=normalize_nt_sequence(
                        seq_to_add,
                        force_upper_case=True,
                        convert_u_into_t=True,
                        convert_non_agct_to_n=True,
                        n_operation="random_assign",
                    ),
                    coordinate=SimpleCoordinate(
                        contig=_transcript_to_use.seqname,
                        strand=strand_repr(_transcript_to_use.strand),
                        start=_transcript_to_use.start,
                        end=_transcript_to_use.end,
                    ),
                )

        def add_transposon(_transposon_to_use: Transcript) -> Optional[SimpleTE]:
            repeat_match_start = int(
                _transposon_to_use.attribute_get("repeat_match_start")
            )
            repeat_match_end = int(_transposon_to_use.attribute_get("repeat_match_end"))
            try:
                seq_to_add = tedb.seq(_transposon_to_use.gene_id)[
                    repeat_match_start:repeat_match_end
                ]
            except KeyError:
                return None

            if len(seq_to_add) < minimal_transposon_len:
                return None
            else:
                return SimpleTE(
                    src_te_name=_transposon_to_use.transcript_id,
                    seq=normalize_nt_sequence(
                        seq_to_add,
                        force_upper_case=True,
                        convert_u_into_t=True,
                        convert_non_agct_to_n=True,
                        n_operation="random_assign",
                    ),
                    coordinate=SimpleCoordinate(
                        contig=_transposon_to_use.seqname,
                        strand=strand_repr(_transposon_to_use.strand),
                        start=_transposon_to_use.start,
                        end=_transposon_to_use.end,
                    ),
                )

        def add_fusion(direction: bool) -> bool:
            """
            :param direction: True = 5 prime (TE before mRNA), False = 3 prime (TE after mRNA)
            :return:
            """
            transcript = rdg.choice(gt.transcript_values)
            transposon_search_direction = direction
            if transcript.strand is False:
                transposon_search_direction = not transposon_search_direction
            if transposon_search_direction:
                possible_transposons = transposon_fi.overlap(
                    (
                        (transcript.seqname, transcript.strand),
                        max(transcript.start - 100000, 0),
                        transcript.start,
                    )
                )
            else:
                possible_transposons = transposon_fi.overlap(
                    (
                        (transcript.seqname, transcript.strand),
                        transcript.end,
                        min(
                            transcript.end + 100000,
                            fav.get_chr_length(transcript.seqname),
                        ),
                    )
                )
            # print(transcript)
            # for i in possible_transposons:
            #     print(transposon_gt.get_transcript(i))
            if not possible_transposons:
                return False
            _add_transposon_result = add_transposon(
                transposon_gt.get_transcript(rdg.choice(possible_transposons))
            )
            _add_transcript_result = add_transcript(transcript)

            if _add_transposon_result is None or _add_transcript_result is None:
                return False
            # print(transcript, possible_transposon, direction)
            if direction:
                new_transcript.l.append(_add_transposon_result)
                new_transcript.l.append(_add_transcript_result)
            else:
                new_transcript.l.append(_add_transcript_result)
                new_transcript.l.append(_add_transposon_result)
            return True

        while len(final_simple_transcripts) < n:
            new_transcript = SimpleTranscript(l=[], d=0)
            state = rdg.choices(choices, weights, k=1)[0]
            if state == FusionTypes.FivePrimeFusion and tedb is not None:
                if not add_fusion(True):
                    continue
            elif state == FusionTypes.ThreePrimeFusion and tedb is not None:
                if not add_fusion(False):
                    continue
            elif state == FusionTypes.GeneGeneFusion:
                add_transcript_result1 = add_transcript(
                    rdg.choice(gt.transcript_values)
                )
                add_transcript_result2 = add_transcript(
                    rdg.choice(gt.transcript_values)
                )
                if (
                    add_transcript_result1 is not None
                    and add_transcript_result2 is not None
                ):
                    new_transcript.l.append(add_transcript_result1)
                    new_transcript.l.append(add_transcript_result2)
                else:
                    continue
            elif state == FusionTypes.GeneOnly:
                add_transcript_result = add_transcript(rdg.choice(gt.transcript_values))
                if add_transcript_result is not None:
                    new_transcript.l.append(add_transcript_result)
                else:
                    continue
            elif state == FusionTypes.TransposonOnly and tedb is not None:
                add_transposon_result = add_transposon(
                    transposon_gt.get_transcript((tedb.draw()[0]))
                )
                if add_transposon_result is not None:
                    new_transcript.l.append(add_transposon_result)
                else:
                    continue
            else:
                continue
            if len(new_transcript.seq) < minimal_seq_len:
                continue
            final_simple_transcripts[f"gtt-{len(final_simple_transcripts)}"] = (
                new_transcript
            )
            pbar.update(1)
        if disable_gmm:
            for v in final_simple_transcripts.values():
                v.d = mu
        else:
            for k, v in depth.simulate_gene_level_depth_gmm(
                gene_names=final_simple_transcripts.keys(),
                mu=mu,
                low_cutoff=low_cutoff,
                high_cutoff_ratio=high_cutoff_ratio,
            ).items():
                final_simple_transcripts[k].d = v
        return cls(final_simple_transcripts)
