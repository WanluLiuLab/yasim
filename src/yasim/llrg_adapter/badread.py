import os
from typing import List, Optional

from yasim.llrg_adapter import BaseLLRGAdapter, LLRG_SHELL_ADAPTER_PATH


class BadreadAdapter(BaseLLRGAdapter):
    model_name: str
    """Filename of pre-defined model."""

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            model_name: str,
            exename: Optional[str] = None,
            **kwargs
    ):
        super().__init__(input_fasta, output_fastq_prefix, depth, exename, **kwargs)
        if self.exename is None:
            self.exename = os.path.join(LLRG_SHELL_ADAPTER_PATH, "badread.sh")
        else:
            self.exename = self.exename
        self.model_name = model_name

    def assemble_cmd(self) -> List[str]:
        cmd = [
            self.exename, "simulate",
            "--reference", self.input_fasta,
            "--quantity", f"{self.depth}x",
            "--error_model", self.model_name,
            "--qscore_model", self.model_name,
            "--start_adapter_seq", "",
            "--end_adapter_seq", ""
        ]
        return cmd

    def move_file_after_finish(self):
        """
        This function is passed since badread pours read into stdout.
        """
        return

    def run(self) -> None:
        self.run_simulator_as_process("badread", stdout_filename=self.output_fastq_prefix + ".fq")
