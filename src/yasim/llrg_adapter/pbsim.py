import glob
import os
from typing import List, Final

from yasim.llrg_adapter import BaseLLRGAdapter, automerge

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
                exename,
                "--prefix", self._tmp_dir,
                "--_depth", str(self._depth),
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_ccs"),
                *other_args,
                self._input_fasta
            ]
        else:
            cmd = [
                exename,
                "--prefix", self._tmp_dir,
                "--_depth", str(self._depth),
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_clr"),
                *other_args,
                self._input_fasta
            ]
    """
    is_ccs: bool
    """
    Whether to simulate CCS or CLR reads.
    """
    _tmp_dir: str
    """Prefix for generated temporary files"""

    _llrg_name: Final[str] = "pbsim"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

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
            depth=depth
        )
        self.is_ccs = is_ccs

        cmd = [
            exename,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--depth", str(self._depth),
        ]
        if self.is_ccs:
            cmd.extend([
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_ccs"),
            ])
        else:
            cmd.extend([
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR, "model_qc_clr"),
            ])
        cmd.extend([
            *other_args,
            self._input_fasta
        ])
        self._cmd = cmd

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

    def _rename_file_after_finish_hook(self):
        automerge(glob.glob(os.path.join(self._tmp_dir, "tmp_????.fastq")), self._output_fastq_prefix + ".fq")

    @property
    def is_pair_end(self) -> bool:
        return False
