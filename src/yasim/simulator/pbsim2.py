import glob
import os
from typing import List

from commonutils import ioctl
from yasim.simulator import Simulator

FILE_DIR = os.path.dirname(__file__)


class SimulatePbsim2(Simulator):
    hmm_model: str
    pbsim2_exename: str

    def assemble_cmd(self) -> List[str]:
        cmd = [
            self.pbsim2_exename,
            "--prefix",
            self.tmp_prefix,
            "--id-prefix",
            self.tmp_prefix,
            "--depth",
            str(self.depth),
            "--hmm_model",
            os.path.join(FILE_DIR, "pbsim2_dist", f"{self.hmm_model}.model"),
            self.input_fasta
        ]
        return cmd

    def move_file_after_finish(self):
        counter = 0
        with ioctl.get_writer(self.output_fastq_prefix + ".fq") as writer:
            for filename in glob.glob(self.tmp_prefix + "_*.fastq"):
                with ioctl.get_reader(filename) as reader:
                    while True:
                        line = reader.readline()
                        if not line:
                            break
                        writer.write(line)
                        counter += 1
                ioctl.rm_rf(filename)
                ioctl.rm_rf(os.path.splitext(filename)[0] + ".maf")
                ioctl.rm_rf(os.path.splitext(filename)[0] + ".ref")

    def run(self) -> None:
        self.run_simulator_as_process("pbsim2")

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            hmm_model: str,
            pbsim2_exename: str = "pbsim2",
            **kwargs):
        super().__init__(input_fasta, output_fastq_prefix, depth, **kwargs)
        self.hmm_model = hmm_model
        self.pbsim2_exename = pbsim2_exename
