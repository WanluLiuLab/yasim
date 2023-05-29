"""
pbsim.py -- Wrapper of PBSIM.
"""

__all__ = (
    "PbsimAdapter",
    "PBSIM_DIST_DIR_PATH",
    "patch_frontend_parser"
)

import argparse
import glob
import os

from labw_utils.typing_importer import List, Final, Mapping, Any
from yasim.helper.frontend import patch_frontend_argument_parser
from yasim.llrg_adapter import BaseProcessBasedLLRGAdapter, automerge

PBSIM_DIST_DIR_PATH = os.path.join(os.path.dirname(__file__), "pbsim_dist")
"""
Where pbsim stores its models
"""


class PbsimAdapter(BaseProcessBasedLLRGAdapter):
    """
    Wrapper of PBSIM.

    Cmdline Spec::

        if self.is_ccs:
            cmd = [
                llrg_executable_path,
                "--prefix", self._tmp_dir,
                "--_depth", str(self._depth),
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR_PATH, "model_qc_ccs"),
                *other_args,
                self._src_fasta_file_path
            ]
        else:
            cmd = [
                llrg_executable_path,
                "--prefix", self._tmp_dir,
                "--_depth", str(self._depth),
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR_PATH, "model_qc_clr"),
                *other_args,
                self._src_fasta_file_path
            ]
    """
    llrg_name: Final[str] = "pbsim"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

    @staticmethod
    def validate_params(**kwargs) -> Mapping[str, Any]:
        return {}

    def __init__(
            self,
            *,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            depth: int,
            llrg_executable_path: str,
            is_trusted: bool,
            is_ccs: bool,
            preserve_intermediate_files: bool,
            other_args: List[str]
    ):
        """
        Initializer.

        :param src_fasta_file_path: Path of source FASTA.
        :param dst_fastq_file_prefix: Prefix od destination FASTQ.
        :param depth: Targeted sequencing depth. Would NOT be related to actual sequencing depth!
        :param llrg_executable_path: Path to LLRG Executable.
        :param other_args: Other arguments to be appended at the bottom of assembled CMD.
        :param is_trusted: Whether to skip input validation test.
        :param is_ccs: Whether to simulate CCS (True) or CLR (False) data.
        :raise LLRGInitializationException: On error.
        """
        super().__init__(
            src_fasta_file_path=src_fasta_file_path,
            dst_fastq_file_prefix=dst_fastq_file_prefix,
            depth=depth,
            llrg_executable_path=llrg_executable_path,
            is_trusted=is_trusted,
            preserve_intermediate_files=preserve_intermediate_files
        )

        cmd = [
            llrg_executable_path,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--depth", str(self._depth),
        ]
        if is_ccs:
            cmd.extend([
                "--data-type", "CCS",
                "--model_qc", os.path.join(PBSIM_DIST_DIR_PATH, "model_qc_ccs"),
            ])
        else:
            cmd.extend([
                "--data-type", "CLR",
                "--model_qc", os.path.join(PBSIM_DIST_DIR_PATH, "model_qc_clr"),
            ])
        cmd.extend([
            *other_args,
            self._src_fasta_file_path
        ])
        self._cmd = cmd

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

    def _post_execution_hook(self):
        automerge(glob.glob(os.path.join(self._tmp_dir, "tmp_????.fastq")), self._dst_fastq_file_prefix + ".fq")

    @property
    def is_pair_end(self) -> bool:
        return False


def patch_frontend_parser(
        parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Patch argument parser with ART arguments.
    """
    parser.add_argument('-c', '--ccs', required=False, help="Simulate CCS instead of CLR", action='store_true')
    parser = patch_frontend_argument_parser(parser, "--preserve_intermediate_files")
    return parser
