import os
import shutil
from typing import List, Union, Final

from labw_utils.commonutils.io import file_system, get_reader, get_writer

from yasim.llrg_adapter import BaseLLRGAdapter


class DwgsimAdapter(BaseLLRGAdapter):
    """
    Wrapper of DWGSIM.

    Cmdline Specs::

        cmd = [
            self.exename,
            "-C", str(self.depth),
            *self.other_args,
            self.input_fasta,
            self.tmp_dir
        ]
    """
    _llrg_name: Final[str] = "dwgsim"
    _require_integer_depth: Final[bool] = False

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
            depth=depth,
            exename=exename,
            other_args=other_args
        )
        self.tmp_dir = self.output_fastq_prefix + ".tmp.d"
        os.makedirs(self.tmp_dir, exist_ok=True)

    def _assemble_cmd_hook(self) -> List[str]:
        cmd = [
            self.exename,
            "-C", str(self.depth),
            *self.other_args,
            self.input_fasta,
            os.path.join(self.tmp_dir, "tmp")
        ]
        return cmd

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
            self.lh.error(f"Unable to find output")

    def run(self) -> None:
        self.run_simulator_as_process()
