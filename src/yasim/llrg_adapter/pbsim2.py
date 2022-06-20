import glob
import os
from typing import List, Optional

from commonutils import shell_utils
from commonutils.io.safe_io import get_writer, get_reader
from yasim.llrg_adapter import BaseLLRGAdapter

PBSIM2_DIST = os.path.join(os.path.dirname(__file__), "pbsim2_dist")
"""
Where PBSIM2 stores its model.
"""


class Pbsim2Adapter(BaseLLRGAdapter):
    hmm_model: str

    def assemble_cmd(self) -> List[str]:
        cmd = [
            self.exename,
            "--prefix", self.tmp_prefix,
            "--id-prefix", self.tmp_prefix,
            "--depth", str(self.depth),
            "--hmm_model", os.path.join(PBSIM2_DIST, f"{self.hmm_model}.model"),
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
                        writer.write(line)
                        counter += 1
                shell_utils.rm_rf(filename)
                shell_utils.rm_rf(os.path.splitext(filename)[0] + ".maf")
                shell_utils.rm_rf(os.path.splitext(filename)[0] + ".ref")

    def run(self) -> None:
        self.run_simulator_as_process("pbsim2")

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            hmm_model: str,
            exename: Optional[str] = None,
            **kwargs):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth,
            exename=exename,
            **kwargs
        )
        self.hmm_model = hmm_model
        if self.exename is None:
            self.exename = "pbsim2"
        else:
            self.exename = exename
