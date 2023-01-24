import os
import shutil
from typing import List, Final

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory

from yasim.llrg_adapter import BaseLLRGAdapter

PBSIM3_DIST = os.path.join(os.path.dirname(__file__), "pbsim3_dist")
"""
Where pbsim3 stores its models
"""


class Pbsim3Adapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM3.

    CMDline Spec::

        cmd = [
            self.exename,
            "--strategy", "trans",
            "--method", "errhmm",
            "--errhmm", self.hmm_model,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--id-prefix", os.path.join(self.tmp_dir, "tmp"),
            "--transcript", self.transcript_tsv_path,
        ]
    """
    hmm_model: str
    """Absolute path to or name of HMM filename"""

    tmp_dir: str
    """Prefix for generated temporary files"""

    _llrg_name: Final[str] = "pbsim3"
    _require_integer_depth: Final[bool] = False

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
        elif os.path.exists(os.path.join(PBSIM3_DIST, f"ERRHMM-{hmm_model}.model")):
            hmm_model = os.path.join(PBSIM3_DIST, f"ERRHMM-{hmm_model}.model")
        else:
            raise ValueError(f"HMM Model {hmm_model} cannot be resolved!")
        self.hmm_model = hmm_model
        self.tmp_dir = self.output_fastq_prefix + ".tmp.d"
        os.makedirs(self.tmp_dir, exist_ok=True)
        self.transcript_tsv_path = os.path.join(
            self.tmp_dir, "transcript.tsv"
        )
        with open(self.transcript_tsv_path, "w") as transcript_writer, \
                FastaViewFactory(
                    filename=input_fasta,
                    read_into_memory=True
                ) as ff:
            transcript_id = ff.chr_names[0]
            transcript_writer.write("\t".join((
                transcript_id,
                str(self.depth),  # Forward
                "0",  # Reverse
                ff.sequence(transcript_id)
            )) + "\n")

    def run(self):
        self.run_simulator_as_process()

    def _assemble_cmd_hook(self) -> List[str]:
        cmd = [
            self.exename,
            "--strategy", "trans",
            "--method", "errhmm",
            "--errhmm", self.hmm_model,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--id-prefix", os.path.join(self.tmp_dir, "tmp"),
            "--transcript", self.transcript_tsv_path,
        ]
        return cmd

    def _rename_file_after_finish_hook(self):
        with open(self.output_fastq_prefix + ".fq", "wb") as writer:
            with open(os.path.join(self.tmp_dir, "tmp.fastq"), "rb") as reader:
                shutil.copyfileobj(reader, writer)
