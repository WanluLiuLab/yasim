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
        if not self.model_name in ("verybad", "verynice"):
            cmd = [
                self.exename, "simulate",
                "--reference", self.input_fasta,
                "--quantity", f"{self.depth}x",
                "--error_model", self.model_name,
                "--qscore_model", self.model_name,
                "--start_adapter_seq", "",
                "--end_adapter_seq", ""
            ]
        elif self.model_name == "verybad":
            cmd = [
                self.exename, "simulate",
                "--reference", self.input_fasta,
                "--quantity", f"{self.depth}x",
                "--glitches", "1000,100,100",
                "--junk_reads", "5",
                "--random_reads", "5",
                "--chimeras", "10",
                "--identity", "75,90,8",
                "--start_adapter_seq", "",
                "--end_adapter_seq", ""
            ]
        elif self.model_name == "verynice":
            cmd = [
                self.exename, "simulate",
                "--reference", self.input_fasta,
                "--quantity", f"{self.depth}x",
                "--error_model", "random",
                "--qscore_model", "ideal",
                "--glitches", "0,0,0",
                "--junk_reads", "0",
                "--random_reads", "0",
                "--chimeras", "0",
                "--identity", "95,100,4",
                "--start_adapter_seq", "",
                "--end_adapter_seq", ""
            ]
        else:
            raise ValueError
        return cmd

    def move_file_after_finish(self):
        """
        This function is passed since badread pours read into stdout.
        """
        return

    def run(self) -> None:
        self.run_simulator_as_process("badread", stdout_filename=self.output_fastq_prefix + ".fq")
