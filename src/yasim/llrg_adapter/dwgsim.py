import os
from typing import List, Union, Final

from labw_utils.commonutils.io import file_system
from yasim.llrg_adapter import BaseLLRGAdapter, LLRGException, autocopy


class DwgsimAdapter(BaseLLRGAdapter):
    """
    Wrapper of DWGSIM.

    Cmdline Specs::

        cmd = [
            exename,
            "-C", str(self._depth),
            *other_args,
            self._input_fasta,
            self._tmp_dir
        ]
    """
    _llrg_name: Final[str] = "dwgsim"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: Union[int, float],
            exename: str,
            other_args: List[str],
    ):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth
        )
        self.tmp_dir = self._output_fastq_prefix + ".tmp.d"

        self._cmd = [
            exename,
            "-C", str(self._depth),
            *other_args,
            self._input_fasta,
            os.path.join(self.tmp_dir, "tmp")
        ]

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

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
            autocopy(r1_file_path, self._output_fastq_prefix + "_1.fq")
            autocopy(r2_file_path, self._output_fastq_prefix + "_2.fq")
            break
        else:
            raise LLRGException("Unable to find output")

    @property
    def is_pair_end(self) -> bool:
        return True
