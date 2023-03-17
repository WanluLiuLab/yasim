import glob
import os
import subprocess
import time
import uuid
from typing import Final, List

import jinja2
import pysam

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.io import file_system
from labw_utils.commonutils.io.safe_io import get_writer
from yasim.llrg_adapter import BaseLLRGAdapter, LLRGException, autocopy, automerge

PBSIM3_DIST = os.path.join(os.path.dirname(__file__), "pbsim3_dist")
"""
Where pbsim3 stores its models
"""

PBSIM3_STRATEGY = ("wgs", "trans")
JINJA2_ENV = jinja2.Environment(loader=jinja2.PackageLoader('yasim.llrg_adapter', 'templates'))
PACB_TEMPLATE = JINJA2_ENV.get_template('pbsim_xml_template.xml')


class Pbsim3Adapter(BaseLLRGAdapter):
    """
    Wrapper of PBSIM3.

    CMDline Spec::

        cmd = [
            exename,
            "--strategy", "trans",
            "--method", self.hmm_method,
            f"--{self.hmm_method}", self.hmm_model,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--transcript" if self.strategy == "trans" else "--genome", self._input_path,
            "--pass-num", str(self.ccs_pass),
            *other_args
        ]
    """
    _tmp_dir: str
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
    _capture_stdout: Final[bool] = False

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
            depth=depth
        )
        self.samtools_path = samtools_path
        self.ccs_path = ccs_path
        self.ccs_pass = ccs_pass
        self.ccs_num_threads = ccs_num_threads

        possible_hmm_model_path = os.path.join(
            PBSIM3_DIST,
            f"{hmm_method.upper()}-{hmm_model}.model"
        )
        if os.path.exists(hmm_model):
            hmm_model = hmm_model
        elif os.path.exists(possible_hmm_model_path):
            hmm_model = possible_hmm_model_path
        else:
            raise ValueError(f"HMM Model {hmm_model} cannot be resolved!")

        if strategy not in PBSIM3_STRATEGY:
            raise ValueError(f"strategy {strategy} should be in {PBSIM3_STRATEGY}!")
        self.strategy = strategy

        if self.strategy == "trans":
            self._input_path = os.path.join(
                self._tmp_dir, "transcript.tsv"
            )
        else:
            self._input_path = self._input_fasta
        self._cmd = [
            exename,
            "--strategy", self.strategy,
            "--method", hmm_method,
            f"--{hmm_method}", hmm_model,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--id-prefix", f"movie{uuid.uuid4()}",
            "--transcript" if self.strategy == "trans" else "--genome", self._input_path,
            "--pass-num", str(self.ccs_pass),
            *other_args
        ]

    def _pre_execution_hook(self) -> None:
        if self.strategy == "trans":
            try:
                with get_writer(self._input_path) as transcript_writer, \
                        FastaViewFactory(
                            filename=self._input_fasta,
                            read_into_memory=True,
                            show_tqdm=False
                        ) as ff:
                    transcript_id = ff.chr_names[0]
                    sequence = ff.sequence(transcript_id)
                    transcript_writer.write("\t".join((
                        transcript_id,
                        str(self._depth),  # Forward
                        "0",  # Reverse
                        sequence
                    )) + "\n")
            except (KeyError, OSError, PermissionError, IndexError) as e:
                raise LLRGException(f"Sequence {transcript_id} from file {self._input_fasta} failed!") from e

    def _rename_file_after_finish_hook(self):
        if self.ccs_pass == 1:
            if self.strategy == "wgs":
                automerge(glob.glob(os.path.join(self._tmp_dir, "tmp_????.fastq")), self._output_fastq_prefix + ".fq")
            else:
                autocopy(os.path.join(self._tmp_dir, "tmp.fastq"), self._output_fastq_prefix + ".fq")
        else:
            subreads_bam_path = os.path.join(self._tmp_dir, "tmp.subreads.bam")
            subreads_sam_path = os.path.join(self._tmp_dir, "tmp.sam")
            subreads_xml_path = os.path.join(self._tmp_dir, "tmp.subreads.xml")
            if not file_system.file_exists(subreads_sam_path):
                return
            with pysam.AlignmentFile(subreads_sam_path, check_sq=False) as sam_reader:
                if len(list(sam_reader.fetch())) == 0:
                    return
            with open(os.path.join(self._tmp_dir, "call_ccs.log"), "wb") as log_writer:
                ccs_xml_path = os.path.join(self._tmp_dir, "tmp.ccs.xml")
                if self._exec_subprocess(
                        [
                            self.samtools_path,
                            "view",
                            subreads_sam_path,
                            "-o", subreads_bam_path
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=log_writer,
                        stderr=log_writer
                ) != 0:
                    return
                with get_writer(subreads_xml_path) as writer:
                    timestamp = time.localtime()
                    writer.write(PACB_TEMPLATE.render(
                        timestamp_file=time.strftime("%y-%m-%dT%H:%M:%S", timestamp),
                        timestamp_simple=time.strftime("%y%m%d_%H%m%S", timestamp),
                        bam_filepath=subreads_bam_path,
                        file_uuid=str(uuid.uuid4())
                    ))
                if self._exec_subprocess(
                        [
                            self.ccs_path,
                            "--report-json", os.path.join(self._tmp_dir, "tmp.ccs.report.json"),
                            "--report-file", os.path.join(self._tmp_dir, "tmp.ccs.report.txt"),
                            "--log-level", "INFO",
                            "--log-file", os.path.join(self._tmp_dir, "tmp.ccs.log"),
                            "--num-threads", str(self.ccs_num_threads),
                            subreads_xml_path,
                            ccs_xml_path
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=log_writer,
                        stderr=log_writer
                ) != 0:
                    return
                with open(self._output_fastq_prefix + ".fq", "wb") as writer:
                    if self._exec_subprocess(
                            [
                                self.samtools_path,
                                "fastq",
                                os.path.join(self._tmp_dir, "tmp.ccs.bam")
                            ],
                            stdin=subprocess.DEVNULL,
                            stdout=writer,
                            stderr=log_writer
                    ) != 0:
                        return

    @property
    def is_pair_end(self) -> bool:
        return False
