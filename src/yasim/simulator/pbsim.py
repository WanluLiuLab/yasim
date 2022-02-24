import glob
import os
from typing import List, Optional

from commonutils import shell_utils
from commonutils.io.safe_io import get_reader, get_writer
from yasim.simulator import Simulator, ADAPTER_SHELL_PATH

FILE_DIR = os.path.dirname(__file__)


class SimulatorPbsim(Simulator):
    is_ccs: bool
    pbsim_exename: str

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            pbsim_exename: Optional[str] = None,
            is_ccs: bool = False,
            **kwargs
    ):
        super().__init__(input_fasta, output_fastq_prefix, depth, **kwargs)
        if pbsim_exename is None:
            self.pbsim_exename = os.path.join(ADAPTER_SHELL_PATH, "pbsim.sh")
        else:
            self.pbsim_exename = pbsim_exename
        self.is_ccs = is_ccs

    def assemble_cmd(self) -> List[str]:
        if self.is_ccs:
            cmd = [
                self.pbsim_exename,
                "--prefix", self.tmp_prefix,
                "--depth", str(self.depth),
                "--data-type", "CCS",
                "--model_qc", os.path.join(FILE_DIR, "pbsim_dist", "model_qc_ccs"),
                self.input_fasta
            ]
        else:
            cmd = [
                self.pbsim_exename,
                "--prefix", self.tmp_prefix,
                "--depth", str(self.depth),
                "--data-type", "CLR",
                "--model_qc", os.path.join(FILE_DIR, "pbsim_dist", "model_qc_clr"),
                self.input_fasta
            ]
        return cmd

    def move_file_after_finish(self):
        counter = 0
        with get_writer(self.output_fastq_prefix + ".fq") as writer:
            for filename in glob.glob(self.tmp_prefix + "_*.fastq"):
                with get_reader(filename) as reader:

                    while True:
                        line = reader.readline()
                        if not line:
                            break
                        if counter % 4 == 0:
                            line = f"@{self.output_fastq_prefix}+{counter}\n"
                        elif counter % 2 == 0:
                            line = "+\n"
                        writer.write(line)
                        counter += 1
                shell_utils.rm_rf(filename)
                shell_utils.rm_rf(os.path.splitext(filename)[0] + ".maf")
                shell_utils.rm_rf(os.path.splitext(filename)[0] + ".ref")

    def run(self) -> None:
        self.run_simulator_as_process("pbsim")
