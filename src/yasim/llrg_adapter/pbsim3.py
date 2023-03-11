import os
import shutil
import subprocess
from typing import List, Final, IO, Union

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.io.safe_io import get_writer
from yasim.llrg_adapter import BaseLLRGAdapter, LLRGException

PBSIM3_DIST = os.path.join(os.path.dirname(__file__), "pbsim3_dist")
"""
Where pbsim3 stores its models
"""

PBSIM3_STRATEGY = ("wgs", "trans")


class Pbsim3Adapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM3.

    CMDline Spec::

        cmd = [
            self.exename,
            "--strategy", "trans",
            "--method", self.hmm_method,
            f"--{self.hmm_method}", self.hmm_model,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--transcript" if self.strategy == "trans" else "--genome", self._input_path,
            "--pass-num", str(self.ccs_pass),
            *self.other_args
        ]
    """
    hmm_model: str
    """Absolute path to or name of HMM filename"""

    hmm_method: str
    """Error-based or quality-score-based model"""

    tmp_dir: str
    """Prefix for generated temporary files"""

    ccs_pass: int
    """Number of CCS passes"""

    samtools_path: str
    """Path to Samtools"""

    ccs_path: str
    """Path to PacBio PBCCS"""

    ccs_num_threads: int
    "Number of threads used when invoking PBCCS"

    strategy: str
    """pbsim3 strategy, can be trans or wgs"""

    _input_path: str
    """Path to input transcript.tsv or genome.fasta"""

    _llrg_name: Final[str] = "pbsim3"
    _require_integer_depth: Final[bool] = False
    _capture_stdout : Final[bool] = False

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            strategy: str,
            hmm_method: str,
            hmm_model: str,
            ccs_pass: int,
            exename: str,
            samtools_path: str,
            ccs_path: str,
            ccs_num_threads: int,
            other_args: List[str]
    ):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth,
            exename=exename,
            other_args=other_args
        )
        self.samtools_path = samtools_path
        self.ccs_path = ccs_path
        self.ccs_pass = ccs_pass
        self.ccs_num_threads = ccs_num_threads

        self.hmm_method = hmm_method
        possible_hmm_model_path = os.path.join(
            PBSIM3_DIST,
            f"{self.hmm_method.upper()}-{hmm_model}.model"
        )
        if os.path.exists(hmm_model):
            hmm_model = hmm_model
        elif os.path.exists(possible_hmm_model_path):
            hmm_model = possible_hmm_model_path
        else:
            raise ValueError(f"HMM Model {hmm_model} cannot be resolved!")
        self.hmm_model = hmm_model
        self.tmp_dir = self.output_fastq_prefix + ".tmp.d"

        if strategy not in PBSIM3_STRATEGY:
            raise ValueError(f"strategy {strategy} should be in {PBSIM3_STRATEGY}!")
        self.strategy = strategy

        if self.strategy == "trans":
            self._input_path = os.path.join(
                self.tmp_dir, "transcript.tsv"
            )
        else:
            self._input_path = self.input_fasta

        self._cmd = [
            self.exename,
            "--strategy", self.strategy,
            "--method", self.hmm_method,
            f"--{self.hmm_method}", self.hmm_model,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--transcript" if self.strategy == "trans" else "--genome", self._input_path,
            "--pass-num", str(self.ccs_pass),
            *self.other_args
        ]

    def _pre_execution_hook(self) -> None:
        try:
            os.makedirs(self.tmp_dir, exist_ok=True)
        except OSError as e:
            raise LLRGException(f"Failed to create temporary directory at {self.tmp_dir}") from e

        if self.strategy == "trans":
            try:
                with get_writer(self._input_path) as transcript_writer, \
                        FastaViewFactory(
                            filename=self.input_fasta,
                            read_into_memory=True,
                            show_tqdm=False
                        ) as ff:
                    transcript_id = ff.chr_names[0]
                    sequence = ff.sequence(transcript_id)
                    transcript_writer.write("\t".join((
                        transcript_id,
                        str(self.depth),  # Forward
                        "0",  # Reverse
                        sequence
                    )) + "\n")
            except (KeyError, OSError, PermissionError) as e:
                raise LLRGException(f"Sequence {transcript_id} from file {self.input_fasta} failed!") from e

    def _rename_file_after_finish_hook(self):
        if self.ccs_pass == 1:
            with open(self.output_fastq_prefix + ".fq", "wb") as writer:
                with open(os.path.join(self.tmp_dir, "tmp.fastq"), "rb") as reader:
                    shutil.copyfileobj(reader, writer)
        else:
            with open(os.path.join(self.tmp_dir, "call_ccs.log"), "wb") as log_writer:
                if self._exec_subprocess(
                        [
                            self.samtools_path,
                            "view",
                            os.path.join(self.tmp_dir, "tmp.sam"),
                            "-o", os.path.join(self.tmp_dir, "tmp.subreads.bam")
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=log_writer,
                        stderr=log_writer
                ) != 0:
                    return
                if self._exec_subprocess(
                        [
                            self.ccs_path,
                            "--report-json", os.path.join(self.tmp_dir, "tmp.ccs.report.json"),
                            "--report-file", os.path.join(self.tmp_dir, "tmp.ccs.report.txt"),
                            "--log-level", "INFO",
                            "--log-file", os.path.join(self.tmp_dir, "tmp.ccs.log"),
                            "--num-threads", str(self.ccs_num_threads),
                            os.path.join(self.tmp_dir, "tmp.subreads.bam"),
                            os.path.join(self.tmp_dir, "tmp.ccs.bam")
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=log_writer,
                        stderr=log_writer
                ) != 0:
                    return
                with open(self.output_fastq_prefix + ".fq", "wb") as writer:
                    if self._exec_subprocess(
                            [
                                self.samtools_path,
                                "fastq",
                                os.path.join(self.tmp_dir, "tmp.ccs.bam")
                            ],
                            stdin=subprocess.DEVNULL,
                            stdout=writer,
                            stderr=log_writer
                    ) != 0:
                        return
