import glob
import shutil
from typing import List, Union, Optional

import commonutils.io.file_system
import commonutils.shutil
from commonutils.compressor import gz_decompress
from yasim.simulator import Simulator


class SimulatorDwgsim(Simulator):
    dwgsim_exename: str

    def __init__(self,
                 input_fasta: str,
                 output_fastq_prefix: str,
                 depth: Union[int, float],
                 dwgsim_exename: Optional[str] = None,
                 **kwargs):
        super().__init__(input_fasta, output_fastq_prefix, depth, **kwargs)
        if dwgsim_exename is None:
            self.dwgsim_exename = "dwgsim"
        else:
            self.dwgsim_exename = dwgsim_exename

    def assemble_cmd(self) -> List[str]:
        cmd = [
            self.dwgsim_exename,
            "-1", "140",
            "-2", "140",
            "-C", str(self.depth),
            self.input_fasta,
            self.tmp_prefix
        ]
        return cmd

    def move_file_after_finish(self):
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
            if not commonutils.io.file_system.file_exists(
                    self.tmp_prefix + suffix_r1) or not commonutils.io.file_system.file_exists(
                self.tmp_prefix + suffix_r2):
                continue
            if suffix_r1.endswith(".gz"):
                gz_decompress(
                    self.tmp_prefix + suffix_r1,
                    self.output_fastq_prefix + f"_1.fq"
                )
                gz_decompress(
                    self.tmp_prefix + suffix_r2,
                    self.output_fastq_prefix + f"_2.fq"
                )
            else:  # Need to archive
                shutil.move(
                    self.tmp_prefix + suffix_r1,
                    self.output_fastq_prefix + f"_1.fq"
                )
                shutil.move(
                    self.tmp_prefix + suffix_r2,
                    self.output_fastq_prefix + f"_2.fq"
                )
            break
        else:
            self.lh.error(f"Unable to find output")
        for filename in glob.glob(self.tmp_prefix + ".bfast.fastq*"):
            commonutils.shutil.rm_rf(filename)

    def run(self) -> None:
        self.run_simulator_as_process("dwgsim")
