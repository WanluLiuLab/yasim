import glob
import os
from typing import List, Optional

from commonutils import ioctl
from yasim.simulator import Simulator, ADAPTER_SHELL_PATH


class SimulatorBadread(Simulator):
    model_name: str
    badread_exename:str

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            model_name: str,
            badread_exename: Optional[str] = None,
            **kwargs
    ):
        super().__init__(input_fasta, output_fastq_prefix, depth, **kwargs)
        if badread_exename is None:
            self.badread_exename = os.path.join(ADAPTER_SHELL_PATH, "badread.sh")
        else:
            self.badread_exename = badread_exename
        self.model_name = model_name

    def assemble_cmd(self) -> List[str]:
        """
        TODO
        """
        cmd = ""
        return cmd

    def move_file_after_finish(self):
        """
        This funcyion is passed due to badread pours read into stdout.
        """
        return

    def run(self) -> None:
        self.run_simulator_as_process("badread", stdout_filename=self.output_fastq_prefix + ".fq")
