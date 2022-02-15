import os
import subprocess
import threading
from abc import abstractmethod
from typing import Dict, Any, Optional

from commonutils import ioctl
from commonutils.logger import get_logger


class Simulator(threading.Thread):
    input_fasta: str
    output_fastq_prefix: str
    """
    Prefix for output FASTQs.
    
    For simulator that generates pair end reads,
    the generated read will named {output_fastq_prefix}_1.fq and {output_fastq_prefix}_2.fq.
    
    For simulators that generates single-end reads,
    the generated read will named {output_fastq_prefix}.fq
    """

    depth: int
    kwargs: Dict[str, Any]
    tmp_prefix: str
    """
    Simulator-based temp directory name.
    """

    def __init__(self,
                 input_fasta: str,
                 output_fastq_prefix: str,
                 depth: int,
                 **kwargs
                 ):
        super(Simulator, self).__init__()
        self.input_fasta = input_fasta
        self.output_fastq_prefix = output_fastq_prefix
        self.depth = depth
        self.kwargs = dict(kwargs)
        self.lh = get_logger(__name__)

    def getopt(self, name: str, default: Optional[Any] = None) -> Any:
        if name in self.kwargs.keys():
            return self.kwargs[name]
        else:
            return default

    @abstractmethod
    def assemble_cmd(self):
        pass

    @abstractmethod
    def move_file_after_finish(self):
        pass

    def run_simulator_as_process(self, simulator_name: str):
        tmp_dir = os.path.join(os.path.dirname(self.output_fastq_prefix), f"{simulator_name}_tmp")
        ioctl.mkdir_p(tmp_dir)
        self.tmp_prefix = os.path.join(tmp_dir, os.path.basename(self.output_fastq_prefix))
        cmd = self.assemble_cmd()
        log_filename = self.output_fastq_prefix + ".log"
        log_handler = open(log_filename, "w")
        self.lh.debug(f"Subprocess {' '.join(cmd)} START")
        p = subprocess.Popen(
            cmd,
            stdin=subprocess.DEVNULL,
            stdout=log_handler,
            stderr=log_handler
        )
        retv = p.wait()
        if retv == 0:
            self.lh.debug(f"Subprocess {' '.join(cmd)} FIN")
            ioctl.rm_rf(log_filename)
            self.move_file_after_finish()
        else:
            self.lh.debug(f"Subprocess {' '.join(cmd)} ERR={retv}")
