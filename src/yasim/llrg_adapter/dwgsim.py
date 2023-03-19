"""
dwgsim.py -- Wrapper of DWGSIM.
"""

__all__ = (
    "DwgsimAdapter",
)

import os
from typing import List, Union, Final, Mapping, Any

from labw_utils.commonutils.io import file_system
from yasim.llrg_adapter import BaseLLRGAdapter, autocopy, NoOutputFileException


class DwgsimAdapter(BaseLLRGAdapter):
    """
    Wrapper of DWGSIM.

    Cmdline Specs::

        cmd = [
            llrg_executable_path,
            "-C", str(self._depth),
            *other_args,
            self._src_fasta_file_path,
            self._tmp_dir
        ]
    """
    llrg_name: Final[str] = "dwgsim"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

    @staticmethod
    def validate_params(**kwargs) -> Mapping[str, Any]:
        return {}

    def __init__(
            self,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            depth: Union[int, float],
            llrg_executable_path: str,
            is_trusted: bool,
            other_args: List[str],
    ):
        """
        Initializer.

        :param src_fasta_file_path: Path of source FASTA.
        :param dst_fastq_file_prefix: Prefix od destination FASTQ.
        :param depth: Targeted sequencing depth. Would NOT be related to actual sequencing depth!
        :param llrg_executable_path: Path to LLRG Executable.
        :param other_args: Other arguments to be appended at the bottom of assembled CMD.
        :param is_trusted: Whether to skip input validation test.
        :raise LLRGInitializationException: On error.
        """
        super().__init__(
            src_fasta_file_path=src_fasta_file_path,
            dst_fastq_file_prefix=dst_fastq_file_prefix,
            depth=depth,
            is_trusted=is_trusted,
            llrg_executable_path=llrg_executable_path
        )
        self._cmd = [
            llrg_executable_path,
            "-C", str(self._depth),
            *other_args,
            self._src_fasta_file_path,
            os.path.join(self._tmp_dir, "tmp")
        ]

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

    def _post_execution_hook(self):
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
            r1_file_path = os.path.join(self._tmp_dir, "tmp" + suffix_r1)
            r2_file_path = os.path.join(self._tmp_dir, "tmp" + suffix_r2)

            if not file_system.file_exists(r1_file_path) or \
                    not file_system.file_exists(r2_file_path):
                continue
            autocopy(r1_file_path, self._dst_fastq_file_prefix + "_1.fq")
            autocopy(r2_file_path, self._dst_fastq_file_prefix + "_2.fq")
            break
        else:
            raise NoOutputFileException("Unable to find output")

    @property
    def is_pair_end(self) -> bool:
        return True
