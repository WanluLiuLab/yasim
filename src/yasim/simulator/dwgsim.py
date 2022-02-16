import glob
import shutil
from typing import List, Union, Optional

from commonutils import ioctl
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
        cmd = [self.dwgsim_exename]
        cmd.append("-C")
        cmd.append(str(self.depth))
        cmd.append(self.input_fasta)
        cmd.append(self.tmp_prefix)
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
            if not ioctl.file_exists(self.tmp_prefix + suffix_r1) or not ioctl.file_exists(self.tmp_prefix + suffix_r2):
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
            ioctl.rm_rf(filename)

    def run(self) -> None:
        self.run_simulator_as_process("dwgsim")
