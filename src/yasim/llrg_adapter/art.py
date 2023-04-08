"""
art.py -- ART adapter for Illumina sequencing.
"""

__all__ = (
    "AVAILABLE_ILLUMINA_ART_SEQUENCER",
    "ArtAdapter",
    "patch_frontend_parser"
)

import argparse
import os
from typing import Dict, List, Tuple, Union, Final, Any, Mapping

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.llrg_adapter import BaseProcessBasedLLRGAdapter, autocopy, LLRGInitializationException

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

_lh = get_logger(__name__)


class ArtAdapter(BaseProcessBasedLLRGAdapter):
    """
    Wrapper of ART.

    Cmdline Specs::

        cmd = [
            llrg_executable_path,
            "--fcov", str(self._depth),
            "--in", self._src_fasta_file_path,
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

    llrg_name: Final[str] = "art"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False
    _is_pair_end: bool

    @staticmethod
    def validate_params(
            sequencer_name: str,
            read_length: int,
            pair_end_fragment_length_mean: int,
            pair_end_fragment_length_std: int,
            is_pair_end: bool,
            **kwargs
    ) -> Mapping[str, Any]:
        _ = kwargs
        retd = {}
        try:
            sequencer_profile = AVAILABLE_ILLUMINA_ART_SEQUENCER[sequencer_name]
        except KeyError as e:
            raise LLRGInitializationException(f"Sequencer Profile for sequencer {sequencer_name} not found!") from e
        if read_length not in sequencer_profile[1]:
            read_length_used = sequencer_profile[1][0]
            _lh.warning("Read length %d not allowed, would use default %d.", read_length, read_length_used)
        else:
            read_length_used = read_length
        retd.update({"read_length": read_length_used})
        if is_pair_end:
            if pair_end_fragment_length_mean * pair_end_fragment_length_std == 0:
                raise LLRGInitializationException(
                    "Please set pair_end_fragment_length_mean and pair_end_fragment_length_std in PE simulation"
                )
        else:
            if pair_end_fragment_length_mean != 0 or pair_end_fragment_length_std != 0:
                _lh.warning(
                    "PE Params pair_end_fragment_length_mean pair_end_fragment_length_std "
                    "should not be used in SE mode, ignored."
                )
        return retd

    def __init__(
            self,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            depth: Union[int, float],
            llrg_executable_path: str,
            is_trusted: bool,
            sequencer_name: str,
            read_length: int,
            pair_end_fragment_length_mean: int,
            pair_end_fragment_length_std: int,
            is_pair_end: bool,
            other_args: List[str]
    ):
        """
        Initializer.

        :param src_fasta_file_path: Path of source FASTA.
        :param dst_fastq_file_prefix: Prefix od destination FASTQ.
        :param depth: Targeted sequencing depth. Would NOT be related to actual sequencing depth!
        :param llrg_executable_path: Path to LLRG Executable.
        :param is_trusted: Whether to skip input validation test.
        :param sequencer_name: name of sequencer to simulate.
        :param read_length: Length of read. Should be defined in corresponding sequencer.
        :param pair_end_fragment_length_mean: [PE Only] The mean size of DNA/RNA fragments for paired-end simulations.
        :param pair_end_fragment_length_std: [PE Only] The standard deviation of DNA/RNA fragment size for
            paired-end simulations.
        :param is_pair_end: Whether the simulation should be pair end.
        :param other_args: Other arguments to be appended at the bottom of assembled CMD.
        :raise LLRGInitializationException: On error.
        """
        super().__init__(
            src_fasta_file_path=src_fasta_file_path,
            dst_fastq_file_prefix=dst_fastq_file_prefix,
            depth=depth,
            llrg_executable_path=llrg_executable_path,
            is_trusted=is_trusted
        )

        if not is_trusted:
            validated_params = ArtAdapter.validate_params(
                sequencer_name=sequencer_name,
                read_length=read_length,
                pair_end_fragment_length_mean=pair_end_fragment_length_mean,
                pair_end_fragment_length_std=pair_end_fragment_length_std,
                is_pair_end=is_pair_end
            )
            read_length = validated_params["read_length"]
        self._is_pair_end = is_pair_end
        self._cmd = [
            self._llrg_executable_path,
            "--fcov", str(self._depth),
            "--in", self._src_fasta_file_path,
            "--samout",
            "--noALN",
            "-l", str(read_length),
            "--out", os.path.join(self._tmp_dir, "tmp"),
            *other_args
        ]
        if self._is_pair_end:
            self._cmd.extend([
                "-p",
                "--sdev", str(pair_end_fragment_length_std),
                "--mflen", str(pair_end_fragment_length_mean)
            ])

    def _post_execution_hook(self):
        if self._is_pair_end:
            autocopy(os.path.join(self._tmp_dir, "tmp1.fq"), self._dst_fastq_file_prefix + "_1.fq")
            autocopy(os.path.join(self._tmp_dir, "tmp2.fq"), self._dst_fastq_file_prefix + "_2.fq")
        else:
            autocopy(os.path.join(self._tmp_dir, "tmp.fq"), self._dst_fastq_file_prefix + ".fq")

    @property
    def is_pair_end(self) -> bool:
        return self._is_pair_end


def patch_frontend_parser(
        parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Patch argument parser with ART arguments.
    """
    parser.add_argument(
        "--sequencer_name",
        required=False,
        help="Name of Illumina Sequencer to Simulate: " + ", ".join((
            f"{name} -- {AVAILABLE_ILLUMINA_ART_SEQUENCER[name][0]}"
            for name in AVAILABLE_ILLUMINA_ART_SEQUENCER.keys()
        )),
        nargs='?',
        choices=AVAILABLE_ILLUMINA_ART_SEQUENCER.keys(),
        type=str,
        action='store',
        default="HS25"
    )
    parser.add_argument(
        "--read_length",
        required=False,
        help="Read length. Sequencer -- Read Length Table: " + ", ".join((
            f"{v[0]} -- {v[1]}"
            for v in AVAILABLE_ILLUMINA_ART_SEQUENCER.values()
        )),
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    parser.add_argument(
        '--pair_end_fragment_length_mean',
        required=False,
        help="[PE Only] The mean size of DNA/RNA fragments for paired-end simulations",
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    parser.add_argument(
        '--pair_end_fragment_length_std',
        required=False,
        help="[PE Only] The standard deviation of DNA/RNA fragment size for paired-end simulations.",
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    parser.add_argument(
        '--is_pair_end',
        required=False,
        help="Whether to use Pair End (PE) Simulation",
        action='store_true'
    )
    return parser
