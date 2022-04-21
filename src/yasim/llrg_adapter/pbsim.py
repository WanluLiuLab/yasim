import glob
import os
from typing import List, Optional

from commonutils import shell_utils
from commonutils.io.safe_io import get_reader, get_writer
from yasim.llrg_adapter import BaseLLRGAdapter, LLRG_SHELL_ADAPTER_PATH

PBSIM_DIST_DIR = os.path.join(os.path.dirname(__file__), "pbsim_dist")
"""
Where pbsim stores its models
"""


class PbsimAdapter(BaseLLRGAdapter):
    is_ccs: bool
    """
    Whether to simulate CCS or CLR reads.
    """

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            exename: Optional[str] = None,
            is_ccs: bool = False,
            **kwargs
    ):
        super().__init__(input_fasta, output_fastq_prefix, depth, exename, **kwargs)
        if self.exename is None:
            self.exename = os.path.join(LLRG_SHELL_ADAPTER_PATH, "pbsim.sh")
        else:
            self.exename = exename
        self.is_ccs = is_ccs
        self._after_cmd_hooks.append(self.hook_del_maf_after_finish)

    def assemble_cmd(self) -> List[str]:
        if self.is_ccs:
            cmd = [
                self.exename,
                "--prefix", self.tmp_prefix,
                "--depth", str(self.depth),
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_ccs"),
                self.input_fasta
            ]
        else:
            cmd = [
                self.exename,
                "--prefix", self.tmp_prefix,
                "--depth", str(self.depth),
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_clr"),
                self.input_fasta
            ]
        return cmd

    def hook_del_maf_after_finish(self):
        for filename in glob.glob(self.tmp_prefix + "_*.fastq"):
            shell_utils.rm_rf(os.path.splitext(filename)[0] + ".maf")
            shell_utils.rm_rf(os.path.splitext(filename)[0] + ".ref")


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


    def run(self) -> None:
        self.run_simulator_as_process("pbsim")
