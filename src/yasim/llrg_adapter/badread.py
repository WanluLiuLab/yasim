from typing import List, Final

from yasim.llrg_adapter import BaseLLRGAdapter


class BadReadAdapter(BaseLLRGAdapter):
    """
    Wrapper for BadRead

    Cmdline Specs::

        if self.model_name == "verybad":
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
                "--end_adapter_seq", "",
                *self.other_args
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
                "--end_adapter_seq", "",
                *self.other_args
            ]
        else:
            cmd = [
                self.exename, "simulate",
                "--reference", self.input_fasta,
                "--quantity", f"{self.depth}x",
                "--error_model", self.model_name,
                "--qscore_model", self.model_name,
                "--start_adapter_seq", "",
                "--end_adapter_seq", "",
                *self.other_args
            ]
    """
    model_name: str
    """Filename of pre-defined model."""

    _llrg_name: Final[str] = "badread"
    _require_integer_depth: Final[bool] = True
    _capture_stdout : Final[bool] = True

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: int,
            model_name: str,
            exename: str,
            other_args: List[str]
    ):
        super().__init__(
            input_fasta=input_fasta,
            output_fastq_prefix=output_fastq_prefix,
            depth=depth,
            exename=exename,
            other_args=other_args
        )

        self.model_name = model_name
        cmd = [
            self.exename, "simulate",
            "--reference", self.input_fasta,
            "--quantity", f"{self.depth}x",
        ]
        if self.model_name == "verybad":
            cmd.extend([
                "--glitches", "1000,100,100",
                "--junk_reads", "5",
                "--random_reads", "5",
                "--chimeras", "10",
                "--identity", "75,90,8"
            ])
        elif self.model_name == "verynice":
            cmd.extend([
                "--error_model", "random",
                "--qscore_model", "ideal",
                "--glitches", "0,0,0",
                "--junk_reads", "0",
                "--random_reads", "0",
                "--chimeras", "0",
                "--identity", "95,100,4"
            ])
        else:
            cmd.extend([
                "--error_model", self.model_name,
                "--qscore_model", self.model_name,
            ])
        cmd.extend([
            "--start_adapter_seq", "",
            "--end_adapter_seq", "",
            *self.other_args
        ])
        self._cmd = cmd

    def _rename_file_after_finish_hook(self):
        """
        This function is passed since badread pours read into stdout.
        """
        pass


    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass
