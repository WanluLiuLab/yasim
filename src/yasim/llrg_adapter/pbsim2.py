import glob
import os
import shutil
from typing import List, Final

from yasim.llrg_adapter import BaseLLRGAdapter, LLRGException

PBSIM2_DIST = os.path.join(os.path.dirname(__file__), "pbsim2_dist")
"""
Where PBSIM2 stores its model.
"""


class Pbsim2Adapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM2.

    CMDline Spec::

        cmd = [
            self.exename,
            "--prefix", self.tmp_dir,
            "--depth", str(self.depth),
            "--hmm_model", self.hmm_model,
            *self.other_args,
            self.input_fasta
        ]
    """
    hmm_model: str
    """Absolute path to or name of HMM filename"""

    tmp_dir: str
    """Prefix for generated temporary files"""

    _llrg_name: Final[str] = "pbsim2"
    _require_integer_depth: Final[bool] = False
    _capture_stdout : Final[bool] = False

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            hmm_model: str,
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
        if os.path.exists(hmm_model):
            hmm_model = hmm_model
        elif os.path.exists(os.path.join(PBSIM2_DIST, f"{hmm_model}.model")):
            hmm_model = os.path.join(PBSIM2_DIST, f"{hmm_model}.model")
        else:
            raise ValueError(f"HMM Model {hmm_model} cannot be resolved!")
        self.hmm_model = hmm_model
        self.tmp_dir = self.output_fastq_prefix + ".tmp.d"
        self._cmd = [
            self.exename,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--depth", str(self.depth),
            "--hmm_model", self.hmm_model,
            *self.other_args,
            self.input_fasta
        ]

    def _pre_execution_hook(self) -> None:
        try:
            os.makedirs(self.tmp_dir, exist_ok=True)
        except OSError as e:
            raise LLRGException(f"Failed to create temporary directory at {self.tmp_dir}") from e

    def _rename_file_after_finish_hook(self):
        with open(self.output_fastq_prefix + ".fq", "wb") as writer:
            for filename in glob.glob(os.path.join(self.tmp_dir, "tmp_????.fastq")):
                with open(filename, "rb") as reader:
                    shutil.copyfileobj(reader, writer)
