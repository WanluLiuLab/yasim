import glob
import os
from typing import List, Final

from yasim.llrg_adapter import BaseLLRGAdapter, automerge

PBSIM2_DIST = os.path.join(os.path.dirname(__file__), "pbsim2_dist")
"""
Where PBSIM2 stores its model.
"""


class Pbsim2Adapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM2.

    CMDline Spec::

        cmd = [
            exename,
            "--prefix", self._tmp_dir,
            "--_depth", str(self._depth),
            "--hmm_model", self.hmm_model,
            *other_args,
            self._input_fasta
        ]
    """
    hmm_model: str
    """Absolute path to or name of HMM filename"""

    _tmp_dir: str
    """Prefix for generated temporary files"""

    _llrg_name: Final[str] = "pbsim2"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

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
            depth=depth
        )
        if os.path.exists(hmm_model):
            hmm_model = hmm_model
        elif os.path.exists(os.path.join(PBSIM2_DIST, f"{hmm_model}.model")):
            hmm_model = os.path.join(PBSIM2_DIST, f"{hmm_model}.model")
        else:
            raise ValueError(f"HMM Model {hmm_model} cannot be resolved!")
        self.hmm_model = hmm_model
        self.tmp_dir = self._output_fastq_prefix + ".tmp.d"
        self._cmd = [
            exename,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--_depth", str(self._depth),
            "--hmm_model", self.hmm_model,
            *other_args,
            self._input_fasta
        ]

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

    def _rename_file_after_finish_hook(self):
        automerge(glob.glob(os.path.join(self.tmp_dir, "tmp_????.fastq")), self._output_fastq_prefix + ".fq")

    @property
    def is_pair_end(self) -> bool:
        return False
