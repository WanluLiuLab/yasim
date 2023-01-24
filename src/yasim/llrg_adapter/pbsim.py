import glob
import os
import shutil
from typing import List, Final

from yasim.llrg_adapter import BaseLLRGAdapter

PBSIM_DIST_DIR = os.path.join(os.path.dirname(__file__), "pbsim_dist")
"""
Where pbsim stores its models
"""


class PbsimAdapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM.

    Cmdline Spec::

        if self.is_ccs:
            cmd = [
                self.exename,
                "--prefix", self.tmp_dir,
                "--depth", str(self.depth),
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_ccs"),
                *self.other_args,
                self.input_fasta
            ]
        else:
            cmd = [
                self.exename,
                "--prefix", self.tmp_dir,
                "--depth", str(self.depth),
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_clr"),
                *self.other_args,
                self.input_fasta
            ]
    """
    is_ccs: bool
    """
    Whether to simulate CCS or CLR reads.
    """
    tmp_dir: str
    """Prefix for generated temporary files"""

    _llrg_name: Final[str] = "pbsim"
    _require_integer_depth: Final[bool] = False

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            is_ccs: bool,
            exename: str,
            other_args: List[str]
    ):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth,
            exename=exename,
            other_args=other_args
        )
        self.is_ccs = is_ccs
        self.tmp_dir = self.output_fastq_prefix + ".tmp.d"
        os.makedirs(self.tmp_dir, exist_ok=True)

    def _assemble_cmd_hook(self) -> List[str]:
        if self.is_ccs:
            cmd = [
                self.exename,
                "--prefix", os.path.join(self.tmp_dir, "tmp"),
                "--depth", str(self.depth),
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_ccs"),
                *self.other_args,
                self.input_fasta
            ]
        else:
            cmd = [
                self.exename,
                "--prefix", os.path.join(self.tmp_dir, "tmp"),
                "--depth", str(self.depth),
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_clr"),
                *self.other_args,
                self.input_fasta
            ]
        return cmd

    def _rename_file_after_finish_hook(self):
        with open(self.output_fastq_prefix + ".fq", "wb") as writer:
            for fn in glob.glob(os.path.join(self.tmp_dir, "tmp_????.fastq")):
                with open(fn, "rb") as reader:
                    shutil.copyfileobj(reader, writer)

    def run(self) -> None:
        self.run_simulator_as_process()
