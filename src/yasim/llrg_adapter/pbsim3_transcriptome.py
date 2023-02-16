import os
import shutil
import subprocess
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
            "--method", self.hmm_method,
            f"--{self.hmm_method}", self.hmm_model,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--transcript", self.transcript_tsv_path,
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

    _llrg_name: Final[str] = "pbsim3"
    _require_integer_depth: Final[bool] = False

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
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
        os.makedirs(self.tmp_dir, exist_ok=True)
        self.transcript_tsv_path = os.path.join(
            self.tmp_dir, "transcript.tsv"
        )
        with open(self.transcript_tsv_path, "w") as transcript_writer, \
                FastaViewFactory(
                    filename=input_fasta,
                    read_into_memory=True,
                    show_tqdm=False
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
            "--method", self.hmm_method,
            f"--{self.hmm_method}", self.hmm_model,
            "--prefix", os.path.join(self.tmp_dir, "tmp"),
            "--transcript", self.transcript_tsv_path,
            "--pass-num", str(self.ccs_pass),
            *self.other_args
        ]
        return cmd

    def _rename_file_after_finish_hook(self):
        if self.ccs_pass == 1:
            with open(self.output_fastq_prefix + ".fq", "wb") as writer:
                with open(os.path.join(self.tmp_dir, "tmp.fastq"), "rb") as reader:
                    shutil.copyfileobj(reader, writer)
        else:
            with open(os.path.join(self.tmp_dir, "call_ccs.log"), "wb") as log_writer:
                subprocess.Popen(
                    [
                        self.samtools_path,
                        "view",
                        os.path.join(self.tmp_dir, "tmp.sam"),
                        "-o", os.path.join(self.tmp_dir, "tmp.subreads.bam")
                    ],
                    stdin=subprocess.DEVNULL,
                    stdout=log_writer,
                    stderr=log_writer
                ).wait()
                subprocess.Popen(
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
                ).wait()
                with open(self.output_fastq_prefix + ".fq", "wb") as writer:
                    subprocess.Popen(
                        [
                            self.samtools_path,
                            "fastq",
                            os.path.join(self.tmp_dir, "tmp.ccs.bam")
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=writer,
                        stderr=log_writer
                    ).wait()
