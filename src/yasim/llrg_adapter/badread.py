"""
badread.py -- Wrapper for BadRead
"""

__all__ = (
    "BadReadAdapter",
    "ALL_POSSIBLE_BADREAD_MODELS",
    "patch_frontend_parser"
)

import argparse

from labw_utils.typing_importer import List, Final, Mapping, Any
from yasim.llrg_adapter import BaseProcessBasedLLRGAdapter, LLRGInitializationException

ALL_POSSIBLE_BADREAD_MODELS = ("nanopore2018", "nanopore2020", "pacbio2016", "verybad", "verynice")
"""All possible badread model names"""


class BadReadAdapter(BaseProcessBasedLLRGAdapter):
    """
    Wrapper for BadRead

    Cmdline Specs::

        if self.model_name == "verybad":
            cmd = [
                llrg_executable_path, "simulate",
                "--reference", self._src_fasta_file_path,
                "--quantity", f"{self._depth}x",
                "--glitches", "1000,100,100",
                "--junk_reads", "5",
                "--random_reads", "5",
                "--chimeras", "10",
                "--identity", "75,90,8",
                "--start_adapter_seq", "",
                "--end_adapter_seq", "",
                *other_args
            ]
        elif self.model_name == "verynice":
            cmd = [
                llrg_executable_path, "simulate",
                "--reference", self._src_fasta_file_path,
                "--quantity", f"{self._depth}x",
                "--error_model", "random",
                "--qscore_model", "ideal",
                "--glitches", "0,0,0",
                "--junk_reads", "0",
                "--random_reads", "0",
                "--chimeras", "0",
                "--identity", "95,100,4",
                "--start_adapter_seq", "",
                "--end_adapter_seq", "",
                *other_args
            ]
        else:
            cmd = [
                llrg_executable_path, "simulate",
                "--reference", self._src_fasta_file_path,
                "--quantity", f"{self._depth}x",
                "--error_model", model_name,
                "--qscore_model", model_name,
                "--start_adapter_seq", "",
                "--end_adapter_seq", "",
                *other_args
            ]
    """

    @staticmethod
    def validate_params(model_name: str, **kwargs) -> Mapping[str, Any]:
        if model_name not in ALL_POSSIBLE_BADREAD_MODELS:
            raise LLRGInitializationException(f"Illegal model name {model_name}")
        return {}

    llrg_name: Final[str] = "badread"
    _require_integer_depth: Final[bool] = True
    _capture_stdout: Final[bool] = True

    def __init__(
            self,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            is_trusted: bool,
            depth: int,
            llrg_executable_path: str,
            model_name: str,
            other_args: List[str]
    ):
        """
        Initializer.

        :param src_fasta_file_path: Path of source FASTA.
        :param dst_fastq_file_prefix: Prefix od destination FASTQ.
        :param depth: Targeted sequencing depth. Would NOT be related to actual sequencing depth!
        :param llrg_executable_path: Path to LLRG Executable.
        :param model_name: Name of pre-defined models. See :py:attr:`ALL_POSSIBLE_BADREAD_MODELS`.
        :param is_trusted: Whether to skip input validation test.
        :param other_args: Other arguments to be appended at the bottom of assembled CMD.
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
            _ = BadReadAdapter.validate_params(model_name=model_name)
        cmd = [
            llrg_executable_path, "simulate",
            "--reference", self._src_fasta_file_path,
            "--quantity", f"{self._depth}x",
        ]
        if model_name == "verybad":
            cmd.extend([
                "--glitches", "1000,100,100",
                "--junk_reads", "5",
                "--random_reads", "5",
                "--chimeras", "10",
                "--identity", "75,90,8"
            ])
        elif model_name == "verynice":
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
                "--error_model", model_name,
                "--qscore_model", model_name,
            ])
        cmd.extend([
            "--start_adapter_seq", "",
            "--end_adapter_seq", "",
            *other_args
        ])
        self._cmd = cmd

    def _post_execution_hook(self):
        """
        This function is passed since badread pours read into stdout.
        """
        pass

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

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
        '--model_name',
        required=True,
        help="Badread model name",
        nargs='?',
        type=str,
        action='store',
        choices=ALL_POSSIBLE_BADREAD_MODELS
    )
    return parser
