import glob
import shutil
from typing import List

from commonutils import ioctl
from yasim.simulator import Simulator


class SimulatorDwgsim(Simulator):
    tmp_prefix: str

    def assemble_cmd(self) -> List[str]:
        cmd = ["dwgsim"]
        cmd.append("-C")
        cmd.append(str(self.depth))
        cmd.append(self.input_fasta)
        cmd.append(self.tmp_prefix)
        return cmd

    def move_file_after_finish(self):
        try:
            shutil.move(
                self.tmp_prefix + ".bwa.read1.fastq.gz",
                self.output_fastq_prefix + "_1.fq.gz"
            )
            shutil.move(
                self.tmp_prefix + ".bwa.read2.fastq.gz",
                self.output_fastq_prefix + "_2.fq.gz"
            )
        except Exception:
            pass
        for filename in glob.glob(self.tmp_prefix + "*"):
            ioctl.rm_rf(filename)

    def run(self) -> None:
        self.run_simulator_as_process("dwgsim")


class SimulatorDwgsimPerfect(SimulatorDwgsim):
    def assemble_cmd(self) -> List[str]:
        cmd = [
            "dwgsim",
            "-d", str(self.getopt("d", 100)),
            "-e", str(self.getopt("e", 0.0)),
            "-E", str(self.getopt("E", 0.0)),
            "-r", str(self.getopt("r", 0.0)),
            "-R", str(self.getopt("R", 0.0)),
            "-F", str(self.getopt("F", 0.0)),
            "-y", str(self.getopt("y", 0.0)),
            "-C", str(self.depth),
            self.input_fasta, self.tmp_prefix
        ]
        return cmd
