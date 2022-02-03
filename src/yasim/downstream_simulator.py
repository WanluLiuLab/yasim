import os
import shutil
import subprocess
import tempfile
import threading
from typing import Dict, Any, Optional, List

from commonutils.logger import get_logger


class Simulator(threading.Thread):
    input_fasta: str
    output_fastq_prefix: str
    depth: int
    kwargs: Dict[str, Any]

    def __init__(self,
                 input_fasta:str,
                 output_fastq_prefix:str,
                 depth:int,
                 **kwargs
                 ):
        super(Simulator, self).__init__()
        self.input_fasta=input_fasta
        self.output_fastq_prefix=output_fastq_prefix
        self.depth=depth
        self.kwargs=dict(kwargs)
        self.lh=get_logger(__name__)

    def getopt(self, name:str, default:Optional[Any]=None) -> Any:
        if name in self.kwargs.keys():
            return self.kwargs[name]
        else:
            return default

class SimulatorDwgsim(Simulator):

    tmp_prefix:str
    def assemble_cmd(self) -> List[str]:
        cmd=["dwgsim"]
        cmd.append("-d")
        cmd.append(str(self.getopt("d", 500)))
        cmd.append("-e")
        cmd.append(str(self.getopt("e", 500)))
        cmd.append("-E")
        cmd.append(str(self.getopt("E", 500)))
        cmd.append("-r")
        cmd.append(str(self.getopt("r", 500)))
        cmd.append("-R")
        cmd.append(str(self.getopt("R", 500)))
        cmd.append("-F")
        cmd.append(str(self.getopt("F", 500)))
        cmd.append("-y")
        cmd.append(str(self.getopt("y", 500)))
        cmd.append("-C")
        cmd.append(str(self.depth))
        cmd.append(self.input_fasta)
        cmd.append(self.tmp_prefix)
        return cmd

    def run(self) -> None:
        self.tmp_prefix = os.path.join(tempfile.mkdtemp(), "F")
        cmd = self.assemble_cmd()
        log_handler=open(self.output_fastq_prefix + ".log", "w")
        self.lh.info(f"Subprocess {' '.join(cmd)} START")
        p = subprocess.Popen(
            cmd,
            stdin=subprocess.DEVNULL,
            stdout=log_handler,
            stderr=log_handler
        )
        retv = p.wait()
        if retv == 0:
            self.lh.info(f"Subprocess {' '.join(cmd)} FIN")
            shutil.move(self.tmp_prefix + ".bwa.read1.fastq", self.output_fastq_prefix + "_1.fq")
            shutil.move(self.tmp_prefix + ".bwa.read2.fastq", self.output_fastq_prefix + "_2.fq")
        else:
            self.lh.info(f"Subprocess {' '.join(cmd)} ERR={retv}")
