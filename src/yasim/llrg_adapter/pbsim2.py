"""
pbsim2.py -- Wrapper of PBSIM2.
"""
__all__ = (
    "Pbsim2Adapter",
    "PBSIM2_DIST_DIR_PATH",
    "PBSIM2_ALL_POSSIBLE_MODELS",
    "patch_frontend_parser"
)

import argparse
import glob
import os
from typing import List, Final, Mapping, Any

from yasim.llrg_adapter import BaseLLRGAdapter, automerge, LLRGInitializationException

PBSIM2_DIST_DIR_PATH = os.path.join(os.path.dirname(__file__), "pbsim2_dist")
"""
Where PBSIM2 stores its model.
"""

PBSIM2_ALL_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename)[0])
    for filename in glob.glob(os.path.join(PBSIM2_DIST_DIR_PATH, "*.model"))
]
"""
Possible PBSIM2 models.
"""


class Pbsim2Adapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM2.

    CMDline Spec::

        cmd = [
            llrg_executable_path,
            "--prefix", self._tmp_dir,
            "--_depth", str(self._depth),
            "--hmm_model", self.hmm_model,
            *other_args,
            self._src_fasta_file_path
        ]
    """
    llrg_name: Final[str] = "pbsim2"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

    @staticmethod
    def validate_params(
            hmm_model: str,
            **kwargs
    ) -> Mapping[str, Any]:
        if os.path.exists(hmm_model):
            hmm_model = hmm_model
        elif os.path.exists(os.path.join(PBSIM2_DIST_DIR_PATH, f"{hmm_model}.model")):
            hmm_model = os.path.join(PBSIM2_DIST_DIR_PATH, f"{hmm_model}.model")
        else:
            raise LLRGInitializationException(f"HMM Model {hmm_model} cannot be resolved!")
        return {"hmm_model": hmm_model}

    def __init__(
            self,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            depth: int,
            llrg_executable_path: str,
            is_trusted: bool,
            hmm_model: str,
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
        :param hmm_model: Name or Path of HMM Model.
        :raise LLRGInitializationException: On error.
        """
        super().__init__(
            src_fasta_file_path=src_fasta_file_path,
            dst_fastq_file_prefix=dst_fastq_file_prefix,
            depth=depth,
            llrg_executable_path=llrg_executable_path,
            is_trusted=is_trusted
        )

        if not is_trusted:
            validated_params = Pbsim2Adapter.validate_params(hmm_model=hmm_model)
            hmm_model = validated_params["hmm_model"]
        self._cmd = [
            llrg_executable_path,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--depth", str(self._depth),
            "--hmm_model", hmm_model,
            *other_args,
            self._src_fasta_file_path
        ]

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
    parser.add_argument(
        '-m',
        '--hmm_model',
        required=True,
        help="Basename or absolute path of HMM file",
        nargs='?',
        type=str,
        action='store',
        choices=PBSIM2_ALL_POSSIBLE_MODELS
    )
    return parser
