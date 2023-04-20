"""
pbsim3.py -- Wrapper of PBSIM3.
"""
__all__ = (
    "Pbsim3Adapter",
    "PBSIM3_DIST_DIR_PATH",
    "PBSIM3_STRATEGY",
    "patch_frontend_parser"
)

import argparse
import glob
import os
import subprocess
import time
import uuid
from labw_utils.typing_importer import Final, List, Mapping, Any, Optional

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.llrg import enhanced_which
from yasim.llrg_adapter import BaseProcessBasedLLRGAdapter, autocopy, automerge, LLRGInitializationException

_lh = get_logger(__name__)

# FIXME: Implement following feature
try:
    import jinja2

    PACB_SUBREAD_XML_TEMPLATE_FILE_PATH = (
        jinja2.
        Environment(
            loader=jinja2.PackageLoader('yasim.llrg_adapter', 'templates'),
            autoescape=True
        ).
        get_template('pbsim_xml_template.xml')
    )
    """PBSIM3 Subread XML Template."""
except ImportError:
    _lh.warning("Jinja2 failed to import. CCS generation will be BAM-based instead of XML-based.")
    jinja2 = None
    PACB_SUBREAD_XML_TEMPLATE_FILE_PATH = None

PBSIM3_DIST_DIR_PATH = os.path.join(os.path.dirname(__file__), "pbsim3_dist")
"""
Where pbsim3 stores its models
"""

PBSIM3_STRATEGY = ("wgs", "trans")
"""PBSIM3 stratergy, can be WGS or Isoform cDNA (trans)"""

PBSIM3_QSHMM_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename.replace("QSHMM-", ""))[0])
    for filename in glob.glob(os.path.join(PBSIM3_DIST_DIR_PATH, "QSHMM-*.model"))
]

PBSIM3_ERRHMM_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename.replace("ERRHMM-", ""))[0])
    for filename in glob.glob(os.path.join(PBSIM3_DIST_DIR_PATH, "ERRHMM-*.model"))
]

_lh = get_logger(__name__)


class Pbsim3Adapter(BaseProcessBasedLLRGAdapter):
    """
    Wrapper of PBSIM3.

    CMDline Spec::

        cmd = [
            llrg_executable_path,
            "--strategy", "trans",
            "--method", self.hmm_method,
            f"--{self.hmm_method}", self.hmm_model,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--transcript" if self.strategy == "trans" else "--genome", self._input_path,
            "--pass-num", str(self._ccs_pass),
            *other_args
        ]
    """
    _ccs_pass: int
    _samtools_executable_path: Optional[str]
    _ccs_executable_path: Optional[str]
    _ccs_num_threads: Optional[int]
    _strategy: str

    _input_path: str
    """Path to input transcript.tsv or genome.fasta"""

    llrg_name: Final[str] = "pbsim3"
    _require_integer_depth: Final[bool] = False
    _capture_stdout: Final[bool] = False

    @staticmethod
    def validate_params(
            hmm_model: str,
            hmm_method: str,
            strategy: str,
            ccs_pass: int,
            ccs_executable_path: Optional[str],
            samtools_executable_path: Optional[str],
            **kwargs
    ) -> Mapping[str, Any]:
        possible_hmm_model_path = os.path.join(
            PBSIM3_DIST_DIR_PATH,
            f"{hmm_method.upper()}-{hmm_model}.model"
        )
        if os.path.exists(hmm_model):
            pass
        elif os.path.exists(possible_hmm_model_path):
            hmm_model = possible_hmm_model_path
        else:
            raise LLRGInitializationException(f"HMM Model {hmm_model} cannot be resolved!")
        if strategy not in PBSIM3_STRATEGY:
            raise LLRGInitializationException(f"strategy {strategy} should be in {PBSIM3_STRATEGY}!")
        if ccs_pass > 1:
            try:
                ccs_executable_path = enhanced_which(ccs_executable_path)
            except FileNotFoundError as e:
                raise LLRGInitializationException(f"PBCCS at {ccs_executable_path} not found!") from e
            try:
                samtools_executable_path = enhanced_which(samtools_executable_path)
            except FileNotFoundError as e:
                raise LLRGInitializationException(f"SAMTOOLS at {samtools_executable_path} not found!") from e
        else:
            if ccs_executable_path is not None:
                _lh.warning("CCS Executable path ignored in CLR mode")
            if samtools_executable_path is not None:
                _lh.warning("SAMTOOLS Executable path ignored in CLR mode")
        return {
            "hmm_model": hmm_model,
            "ccs_executable_path": ccs_executable_path,
            "samtools_executable_path": samtools_executable_path
        }

    def __init__(
            self,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            depth: int,
            llrg_executable_path: str,
            is_trusted: bool,
            strategy: str,
            hmm_method: str,
            hmm_model: str,
            samtools_executable_path: Optional[str],
            ccs_executable_path: Optional[str],
            ccs_num_threads: Optional[int],
            ccs_pass: int,
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
        :param ccs_num_threads: [CCS] Number of threads used when invoking PBCCS.
        :param ccs_executable_path: [CCS] Path to PacBio PBCCS
        :param samtools_executable_path: [CCS] Path to Samtools.
        :param ccs_pass: Number of CCS passes, 1 for CLR, >1 for CCS.
        :param strategy: PBSIM3 stratergy, can be WGS or Isoform cDNA (trans).
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
            validated_params = Pbsim3Adapter.validate_params(
                hmm_model=hmm_model,
                hmm_method=hmm_method,
                strategy=strategy,
                ccs_pass=ccs_pass,
                samtools_executable_path=samtools_executable_path,
                ccs_executable_path=ccs_executable_path
            )
            hmm_model = validated_params["hmm_model"]
        self._strategy = strategy
        self._samtools_executable_path = samtools_executable_path
        self._ccs_executable_path = ccs_executable_path
        self._ccs_pass = ccs_pass
        self._ccs_num_threads = ccs_num_threads

        possible_hmm_model_path = os.path.join(
            PBSIM3_DIST_DIR_PATH,
            f"{hmm_method.upper()}-{hmm_model}.model"
        )
        if os.path.exists(hmm_model):
            pass
        elif os.path.exists(possible_hmm_model_path):
            hmm_model = possible_hmm_model_path
        else:
            raise LLRGInitializationException(f"HMM Model {hmm_model} cannot be resolved!")

        if strategy not in PBSIM3_STRATEGY:
            raise LLRGInitializationException(f"strategy {strategy} should be in {PBSIM3_STRATEGY}!")
        self._strategy = strategy

        if self._strategy == "trans":
            self._input_path = os.path.join(
                self._tmp_dir, "transcript.tsv"
            )
            strategy_params = [
                "--transcript", self._input_path
            ]
        else:
            self._input_path = self._src_fasta_file_path
            strategy_params = [
                "--depth", str(depth),
                "--genome", self._input_path
            ]
        self._cmd = [
            llrg_executable_path,
            "--strategy", self._strategy,
            "--method", hmm_method,
            f"--{hmm_method}", hmm_model,
            "--prefix", os.path.join(self._tmp_dir, "tmp"),
            "--id-prefix", f"movie{uuid.uuid4()}",
            *strategy_params,
            "--pass-num", str(self._ccs_pass),
            *other_args
        ]

    def _pre_execution_hook(self) -> None:
        if self._strategy == "trans":
            try:
                with get_writer(self._input_path) as transcript_writer, \
                        FastaViewFactory(
                            filename=self._src_fasta_file_path,
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
            except (KeyError, OSError, IndexError) as e:
                raise LLRGInitializationException(
                    f"Sequence {transcript_id} from file {self._src_fasta_file_path} failed!") from e

    def _ccs_to_fastq(self, prefix: str):
        subreads_bam_path = os.path.join(self._tmp_dir, f"{prefix}.subreads.bam")
        subreads_sam_path = os.path.join(self._tmp_dir, f"{prefix}.sam")
        subreads_xml_path = os.path.join(self._tmp_dir, f"{prefix}.subreads.xml")
        ccs_bam_path = os.path.join(self._tmp_dir, f"{prefix}.ccs.bam")
        ccs_xml_path = os.path.join(self._tmp_dir, f"{prefix}.ccs.xml")
        output_fastq_path = os.path.join(self._tmp_dir, f"{prefix}.ccs.fq")
        with get_writer(os.path.join(self._tmp_dir, "call_ccs.log"), is_binary=True) as log_writer:
            if self._exec_subprocess(
                    [
                        self._samtools_executable_path,
                        "view",
                        subreads_sam_path,
                        "-o", subreads_bam_path
                    ],
                    stdin=subprocess.DEVNULL,
                    stdout=log_writer,
                    stderr=log_writer
            ) != 0:
                return
            if jinja2 is not None:
                with get_writer(subreads_xml_path) as writer:
                    timestamp = time.localtime()
                    writer.write(PACB_SUBREAD_XML_TEMPLATE_FILE_PATH.render(
                        timestamp_file=time.strftime("%y-%m-%dT%H:%M:%S", timestamp),
                        timestamp_simple=time.strftime("%y%m%d_%H%m%S", timestamp),
                        bam_filepath=subreads_bam_path,
                        file_uuid=str(uuid.uuid4())
                    ))
            if self._exec_subprocess(
                    [
                        self._ccs_executable_path,
                        "--report-json", os.path.join(self._tmp_dir, f"{prefix}.ccs.report.json"),
                        "--report-file", os.path.join(self._tmp_dir, f"{prefix}.ccs.report.txt"),
                        "--log-level", "INFO",
                        "--log-file", os.path.join(self._tmp_dir, f"{prefix}.ccs.log"),
                        "--num-threads", str(self._ccs_num_threads),
                        subreads_xml_path if jinja2 is not None else subreads_bam_path,
                        ccs_xml_path
                    ],
                    stdin=subprocess.DEVNULL,
                    stdout=log_writer,
                    stderr=log_writer
            ) != 0:
                return
            with get_writer(output_fastq_path, is_binary=True) as writer:
                if self._exec_subprocess(
                        [
                            self._samtools_executable_path,
                            "fastq",
                            ccs_bam_path
                        ],
                        stdin=subprocess.DEVNULL,
                        stdout=writer,
                        stderr=log_writer
                ) != 0:
                    return

    def _post_execution_hook(self):
        if self._ccs_pass == 1:
            if self._strategy == "wgs":
                automerge(glob.glob(os.path.join(self._tmp_dir, "tmp_????.fastq")), self._dst_fastq_file_prefix + ".fq")
            else:
                autocopy(os.path.join(self._tmp_dir, "tmp.fastq"), self._dst_fastq_file_prefix + ".fq")
        else:
            if self._strategy == "wgs":
                for fp in glob.glob(os.path.join(self._tmp_dir, "tmp_????.sam")):
                    prefix = os.path.basename(fp).split(".")[0]
                    self._ccs_to_fastq(prefix=prefix)
                automerge(glob.glob(os.path.join(self._tmp_dir, "tmp_????.ccs.fq")),
                          self._dst_fastq_file_prefix + ".fq")
            else:
                self._ccs_to_fastq(prefix="tmp")
                autocopy(os.path.join(self._tmp_dir, "tmp.ccs.fq"), self._dst_fastq_file_prefix + ".fq")

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
        help="Basename of HMM file. "
             f"If you select errhmm in hmm_method, it would be {PBSIM3_ERRHMM_POSSIBLE_MODELS}"
             f"If you select qshmm in hmm_method, it would be {PBSIM3_QSHMM_POSSIBLE_MODELS}",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        "-M",
        "--hmm_method",
        required=True,
        help="Whether to simulate using quality score (as PBSIM2) or error profile (new)",
        nargs='?',
        type=str,
        action='store',
        choices=["errhmm", "qshmm"]
    )
    parser.add_argument(
        "--ccs_pass",
        required=False,
        help="CCS Multipass Settings. Use 1 for CLR and others for CCS.",
        nargs='?',
        type=int,
        action='store',
        default=1
    )
    parser.add_argument(
        '--ccs_path',
        required=False,
        help="Executable name of ccs or pbccs. Omitted if ccs_pass == 1.",
        nargs='?',
        type=str,
        action='store',
        default="ccs"
    )
    parser.add_argument(
        '--samtools_path',
        required=False,
        help="Executable name of samtools. Omitted if ccs_pass == 1.",
        nargs='?',
        type=str,
        action='store',
        default="samtools"
    )
    parser.add_argument(
        '--strategy',
        required=False,
        help="Whether to use transcript (trans) mode or wgs (wgs) mode",
        choices=PBSIM3_STRATEGY,
        type=str,
        action='store',
        default="wgs"
    )
    return parser
