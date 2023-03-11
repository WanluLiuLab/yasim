import os
import shutil
from typing import Dict, List, Tuple, Union, Final

from labw_utils.commonutils.io import file_system, get_reader, get_writer

from yasim.llrg_adapter import BaseLLRGAdapter, LLRGException

AVAILABLE_ILLUMINA_ART_SEQUENCER: Dict[str, Tuple[str, List[int]]] = {
        "GA1":("GenomeAnalyzer I", [36,44]),
        "GA2":("GenomeAnalyzer II", [50, 75]),
        "HS10":("HiSeq 1000", [100]),
        "HS20":("HiSeq 2000", [100]),
        "HS25":("HiSeq 2500", [125, 150]),
        "HSXn":("HiSeqX PCR free", [150]),
        "HSXt":("HiSeqX TruSeq", [150]),
        "MinS":("MiniSeq TruSeq", [50]),
        "MSv1":("MiSeq v1", [250]),
        "MSv3":("MSv3 - MiSeq v3", [250]),
        "NS50":("NextSeq500 v2", [75])
}

class ArtAdapter(BaseLLRGAdapter):
    """
    Wrapper of ART.

    Cmdline Specs::

        cmd = [
            self.exename,
            "-C", str(self.depth),
            *self.other_args,
            self.input_fasta,
            self.tmp_dir
        ]
    """
    _llrg_name: Final[str] = "art"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: Union[int, float],
            exename: str,
            sequencer: str,
            rlen: int,
            is_pair_end: bool,
            other_args: List[str],
    ):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth,
            exename=exename,
            other_args=other_args
        )
        self.tmp_dir = self.output_fastq_prefix + ".tmp.d"
        sequencer_profile = AVAILABLE_ILLUMINA_ART_SEQUENCER[sequencer]

        if rlen not in sequencer_profile[2]:
            rlen = sequencer_profile[2][0]

        self._cmd = [
            self.exename,
            "-C", str(self.depth),
            *self.other_args,
            self.input_fasta,
            
            "-l", rlen,
            "-O", os.path.join(self.tmp_dir, "tmp")
        ]

    def _pre_execution_hook(self) -> None:
        try:
            os.makedirs(self.tmp_dir, exist_ok=True)
        except OSError as e:
            raise LLRGException(f"Failed to create temporary directory at {self.tmp_dir}") from e


    def _rename_file_after_finish_hook(self):
        try_read1_suffix = (
            ".bwa.read1.fastq.gz",
            ".bwa.read1.fastq",
            ".bwa.read1.fq.gz",
            ".bwa.read1.fq"
        )
        try_read2_suffix = (
            ".bwa.read2.fastq.gz",
            ".bwa.read2.fastq",
            ".bwa.read2.fq.gz",
            ".bwa.read2.fq"
        )
        for suffix_r1, suffix_r2 in zip(try_read1_suffix, try_read2_suffix):
            r1_file_path = os.path.join(self.tmp_dir, "tmp" + suffix_r1)
            r2_file_path = os.path.join(self.tmp_dir, "tmp" + suffix_r2)

            if not file_system.file_exists(r1_file_path) or \
                    not file_system.file_exists(r2_file_path):
                continue
            with get_reader(r1_file_path, is_binary=True) as r1, \
                    get_writer(self.output_fastq_prefix + "_1.fq", is_binary=True) as w1:
                shutil.copyfileobj(r1, w1)
            with get_reader(r2_file_path, is_binary=True) as r2, \
                    get_writer(self.output_fastq_prefix + "_2.fq", is_binary=True) as w2:
                shutil.copyfileobj(r2, w2)
            break
        else:
            self._lh.error(f"Unable to find output")
