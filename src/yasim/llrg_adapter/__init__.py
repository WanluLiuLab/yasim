import logging
import subprocess
import threading
from abc import abstractmethod
from typing import Union, List, IO, Optional

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger


class LLRGException(RuntimeError):
    def __init__(self, contents: str):
        self.contents = contents

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return self.contents


class BaseLLRGAdapter(threading.Thread):
    """
    Base class of LLRG Python adapter.

    The LLRG adapter is a thread-safe Python wrapper of various LLRGs.
    They can be used to simulate DNA and RNA data with pre-defined arguments.


    It performs following operations:

    - Assemble the command that calls LLRG.
    - Execute the command using :py:mod:`subprocess`.
    - Move generated files to ``output_fastq_prefix``.
    """

    _lh: logging.Logger
    """
    Class logger handler
    """

    _cmd: Optional[List[str]]
    """Assembled Commandline"""

    input_fasta: str
    """
    Input reference FASTA,
    can be DNA or transcript cDNA,
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

    tmp_dir: str
    """
    Simulator-based temp directory name.
    """

    exename: str
    """
    name of executable.
    Can be name, relative or absolute path to LLRG Shell Adapters or LLRG itself.
    """

    # Following fields are left for LLRGs.
    _llrg_name: str
    """Class attribute, indicating name"""

    _require_integer_depth: bool
    """Class attribute, indicating whether the depts should be converted to ine"""

    _capture_stdout: bool
    """Whether this simulator pours data into stdout"""

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: Union[int, float],
            exename: str,
            other_args: List[str]
    ):
        super(BaseLLRGAdapter, self).__init__()
        if not hasattr(self, "_llrg_name"):
            raise TypeError
        if not hasattr(self, "_require_integer_depth"):
            raise TypeError
        if not hasattr(self, "_capture_stdout"):
            raise TypeError
        self.input_fasta = input_fasta
        self.output_fastq_prefix = output_fastq_prefix
        self.depth = int(depth) if self._require_integer_depth else depth
        self.exename = exename
        self._lh = get_logger(__name__)
        self.other_args = other_args
        self._cmd = None

    @abstractmethod
    def _pre_execution_hook(self) -> None:
        """
        Additional steps to do before starting the simulator process
        """
        raise NotImplementedError

    @abstractmethod
    def _rename_file_after_finish_hook(self):
        """
        Move the file into desired destination.
        May involve (de-)compression.
        """
        raise NotImplementedError

    def run(self):
        if self._cmd is None:
            self._lh.error("Commandline Assembly Failed!")
            return
        try:
            self._pre_execution_hook()
        except LLRGException as e:
            self._lh.error("Exception %s caught at pre-execution time", str(e))
            return

        if self._capture_stdout:
            with open(self.output_fastq_prefix + ".log", "wb") as subprocess_log_handler, \
                    open(self.output_fastq_prefix + ".fq", "wb") as stdout_handler:
                retv = self._exec_subprocess(
                    self._cmd,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_handler,
                    stderr=subprocess_log_handler
                )
        else:
            with open(self.output_fastq_prefix + ".log", "wb") as subprocess_log_handler:
                retv = self._exec_subprocess(
                    self._cmd,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess_log_handler,
                    stderr=subprocess_log_handler
                )
        if retv != 0:
            return
        try:
            self._rename_file_after_finish_hook()
        except LLRGException as e:
            self._lh.error("Exception %s caught at file renaming", str(e))
            return

    def _exec_subprocess(
            self,
            cmd: List[str],
            stdin: Union[IO, int],
            stdout: Union[IO, int],
            stderr: Union[IO, int]
    ) -> int:
        """Wrapper of :py:class:`subprocess.Popen` which logs."""
        self._lh.debug(f"Subprocess {' '.join(cmd)} START")
        retv = subprocess.Popen(
            cmd,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr
        ).wait()
        if retv == 0:
            self._lh.debug(f"Subprocess {' '.join(cmd)} FIN")
        else:
            self._lh.error(f"Subprocess {' '.join(cmd)} ERR={retv}")
        return retv
