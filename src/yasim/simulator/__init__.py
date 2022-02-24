import os
import subprocess
import threading
from abc import abstractmethod
from typing import Dict, Any, Optional, Union, List

from commonutils import shell_utils
from commonutils.stdlib_helper.logger_helper import get_logger

ADAPTER_SHELL_PATH = os.path.join(os.path.dirname(__file__), "shell")


class Simulator(threading.Thread):
    input_fasta: str
    """
    Input reference FASTA, can be DNA or transcriptom cDNA.
    """
    output_fastq_prefix: str
    """
    Prefix for output FASTQs.
    
    For simulator that generates pair end reads,
    the generated read will named {output_fastq_prefix}_1.fq and {output_fastq_prefix}_2.fq.
    
    For simulators that generates single-end reads,
    the generated read will named {output_fastq_prefix}.fq
    """

    depth: Union[int, float]
    """
    Sequencing depth or coverage.
    """

    kwargs: Dict[str, Any]
    """
    Other arguments for simulator.
    """

    tmp_prefix: str
    """
    Simulator-based temp directory name.
    """

    def __init__(self,
                 input_fasta: str,
                 output_fastq_prefix: str,
                 depth: Union[int, float],
                 **kwargs
                 ):
        super(Simulator, self).__init__()
        self.input_fasta = input_fasta
        self.output_fastq_prefix = output_fastq_prefix
        self.depth = depth
        self.kwargs = dict(kwargs)
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
        Move the file into desired destination. May involve (de-)compression.
        """
        pass

    def run_simulator_as_process(self, simulator_name: str, stdout_filename: Optional[str] = None) -> int:
        """
        Use a subprocess to run simulator.

        :param simulator_name: Name to use in logs.
        :param stdout_filename: If the simulator creates output in stdout, set this param to destination filename.
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
