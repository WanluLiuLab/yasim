import os
from typing import Dict, List, Tuple, Union, Final

from yasim.llrg_adapter import BaseLLRGAdapter, autocopy, LLRGInitializationException

AVAILABLE_ILLUMINA_ART_SEQUENCER: Dict[str, Tuple[str, List[int]]] = {
    "GA1": ("GenomeAnalyzer I", [36, 44]),
    "GA2": ("GenomeAnalyzer II", [50, 75]),
    "HS10": ("HiSeq 1000", [100]),
    "HS20": ("HiSeq 2000", [100]),
    "HS25": ("HiSeq 2500", [125, 150]),
    "HSXn": ("HiSeqX PCR free", [150]),
    "HSXt": ("HiSeqX TruSeq", [150]),
    "MinS": ("MiniSeq TruSeq", [50]),
    "MSv1": ("MiSeq v1", [250]),
    "MSv3": ("MSv3 - MiSeq v3", [250]),
    "NS50": ("NextSeq500 v2", [75])
}


class ArtAdapter(BaseLLRGAdapter):
    """
    Wrapper of ART.

    Cmdline Specs::

        cmd = [
            exename,
            "--fcov", str(self._depth),
            "--in", self._input_fasta,
            "--samout",
            "--noALN",
            "-l", rlen,
            "-O", os.path.join(self._tmp_dir, "tmp"),
            *other_args
        ]
    """

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

    _llrg_name: Final[str] = "art"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False
    _is_pair_end: bool

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: Union[int, float],
            exename: str,
            sequencer: str,
            rlen: int,
            mflen_mean: int,
            mflen_std: int,
            is_pair_end: bool,
            other_args: List[str],
    ):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth
        )
        self._is_pair_end = is_pair_end
        self._tmp_dir = self._output_fastq_prefix + ".tmp.d"
        try:
            sequencer_profile = AVAILABLE_ILLUMINA_ART_SEQUENCER[sequencer]
        except KeyError as e:
            raise LLRGInitializationException from e

        if rlen not in sequencer_profile[1]:
            rlen = sequencer_profile[1][0]
        self._cmd = [
            exename,
            "--fcov", str(self._depth),
            "--in", self._input_fasta,
            "--samout",
            "--noALN",
            "-l", str(rlen),
            "--out", os.path.join(self._tmp_dir, "tmp"),
            *other_args
        ]
        if self._is_pair_end:
            self._cmd.extend([
                "-p",
                "--sdev", str(mflen_std),
                "--mflen", str(mflen_mean)
            ])

    def _post_execution_hook(self):
        if self._is_pair_end:
            autocopy(os.path.join(self._tmp_dir, "tmp1.fq"), self._output_fastq_prefix + "_1.fq")
            autocopy(os.path.join(self._tmp_dir, "tmp2.fq"), self._output_fastq_prefix + "_2.fq")
        else:
            autocopy(os.path.join(self._tmp_dir, "tmp.fq"), self._output_fastq_prefix + ".fq")

    @property
    def is_pair_end(self) -> bool:
        return self._is_pair_end
