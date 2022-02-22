import os
from typing import List, Optional

from yasim.simulator import Simulator, ADAPTER_SHELL_PATH


class SimulatorSimlord(Simulator):
    simlord_exename: str

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            simlord_exename: Optional[str] = None,
            **kwargs
    ):
        super().__init__(input_fasta, output_fastq_prefix, depth, **kwargs)
        if simlord_exename is None:
            self.simlord_exename = os.path.join(ADAPTER_SHELL_PATH, "simlord.sh")
        else:
            self.simlord_exename = simlord_exename

    def assemble_cmd(self) -> List[str]:
        cmd = [
            self.simlord_exename,
            "--read-reference", self.input_fasta,
            "--coverage", str(self.depth),
            self.output_fastq_prefix
        ]
        return cmd

    def move_file_after_finish(self):
        """
        This funcyion is passed due to simlord pours read into :py:func:`self.output_fastq_prefix`.
        """
        return

    def run(self) -> None:
        self.run_simulator_as_process("simlord")
