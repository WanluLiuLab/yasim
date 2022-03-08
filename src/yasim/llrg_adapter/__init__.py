import os
import subprocess
import threading
from abc import abstractmethod
from typing import Dict, Any, Optional, Union, List

from commonutils import shell_utils
from commonutils.stdlib_helper.logger_helper import get_logger

LLRG_SHELL_ADAPTER_PATH = os.path.join(os.path.dirname(__file__), "shell")


class BaseLLRGAdapter(threading.Thread):
    """
    Base class of LLRG Python adapter.
    See :doc:`/src/devel/design` for more details.

    It performs following operations:

    - Assemble the command that calls LLRG.
    - Execute the command using :py:mod:`subprocess`.
    - Move generated files to ``output_fastq_prefix``.
    """

    input_fasta: str
    """
    Input reference FASTA,
    can be DNA or transcriptom cDNA,
    can contain more than 1 entries.
    """

    output_fastq_prefix: str
    """
    Prefix for output FASTQs.
    
    For LLRG that generates pair end (PE) reads,
    the generated read will named {output_fastq_prefix}_1.fq and {output_fastq_prefix}_2.fq.
    
    For LLRG that generates single-end (SE) reads,
    the generated read will named {output_fastq_prefix}.fq
    """

    depth: Union[int, float]
    """
    Sequencing depth or coverage.
    
    .. warning:: This is NOT final read count!
    """

    kwargs: Dict[str, Any]
    """
    Other arguments for LLRG.
    """

    tmp_prefix: str
    """
    Simulator-based temp directory name.
    """

    exename: str
    """
    name of executable.
    Can be name, relative or absolute path to LLRG Shell Adapters or LLRG itself.
    """

    def __init__(self,
                 input_fasta: str,
                 output_fastq_prefix: str,
                 depth: Union[int, float],
                 exename: Optional[str] = None,
                 **kwargs
                 ):
        super(BaseLLRGAdapter, self).__init__()
        self.input_fasta = input_fasta
        self.output_fastq_prefix = output_fastq_prefix
        self.depth = depth
        self.kwargs = dict(kwargs)
        self.exename = exename
        self.lh = get_logger(__name__)

    @abstractmethod
    def assemble_cmd(self) -> List[str]:
        """
        Assemble the command line into a list of string.
        """
        pass

    @abstractmethod
    def move_file_after_finish(self):
        """
        Move the file into desired destination.
        May involve (de-)compression.
        """
        pass

    def run_simulator_as_process(self, simulator_name: str, stdout_filename: Optional[str] = None) -> int:
        """
        Use a subprocess to run LLRG.

        :param simulator_name: Name to use in logs.
        :param stdout_filename: If the LLRG creates output in stdout, set this param to destination filename.
        :return: The return value of top-level LLRG process. 0 for normal.
        """
        tmp_dir = os.path.join(os.path.dirname(self.output_fastq_prefix), f"{simulator_name}_tmp")
        shell_utils.mkdir_p(tmp_dir)
        self.tmp_prefix = os.path.join(tmp_dir, os.path.basename(self.output_fastq_prefix))
        cmd = self.assemble_cmd()
        log_filename = self.output_fastq_prefix + ".log"
        log_handler = open(log_filename, "w")
        self.lh.debug(f"Subprocess {' '.join(cmd)} START")
        if stdout_filename is None:
            stdout = log_handler
        else:
            stdout = open(stdout_filename, "wb")
        p = subprocess.Popen(
            cmd,
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=log_handler
        )
        retv = p.wait()
        if retv == 0:
            self.lh.debug(f"Subprocess {' '.join(cmd)} FIN")
            shell_utils.rm_rf(log_filename)
            self.move_file_after_finish()
        else:
            self.lh.debug(f"Subprocess {' '.join(cmd)} ERR={retv}")
        return retv
