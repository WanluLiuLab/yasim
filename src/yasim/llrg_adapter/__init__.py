import logging
import os.path
import shutil
import subprocess
import threading
from abc import abstractmethod
from typing import Union, List, IO, Optional, Iterable

from labw_utils.commonutils.io import get_reader, get_writer, file_system
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger


class LLRGException(RuntimeError):
    """Some error raised due to (improperly configured) LLRG"""

    def __init__(self, contents: str):
        self.contents = contents

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return self.contents


def autocopy(in_fn: str, out_fn: str) -> None:
    """
    Copy one file to another place, with automatic extraction of GZipped files

    :param in_fn: Input filename.
    :param out_fn: Output filename.
    :raises LLRGException: If error occurs in copying.
    """
    try:
        with get_reader(in_fn, is_binary=True) as r1, \
                get_writer(out_fn, is_binary=True) as w1:
            shutil.copyfileobj(r1, w1)
    except (FileNotFoundError, OSError, PermissionError, IOError) as e:
        raise LLRGException(f"Copy file {in_fn} -> {out_fn} failed!") from e


def automerge(in_fns: Iterable[str], out_fn: str) -> None:
    """
    Merge multiple files into one. See :py:func:`autocopy`.
    """
    with get_writer(out_fn, is_binary=True) as w1:
        for in_fn in in_fns:
            try:
                with get_reader(in_fn, is_binary=True) as r1:
                    shutil.copyfileobj(r1, w1)
            except (FileNotFoundError, OSError, PermissionError,) as e:
                raise LLRGException(f"Copy file {in_fn} -> {out_fn} failed!") from e


class BaseLLRGAdapter(threading.Thread):
    """
    Base class of LLRG Python adapter.

    The LLRG adapter is a thread-safe Python wrapper of various LLRGs.
    They can be used to simulate DNA and RNA data with pre-defined arguments.


    It performs following operations:

    - Assemble the command that calls LLRG.
    - Execute the command using :py:mod:`subprocess`.
    - Move generated files to ``_output_fastq_prefix``.
    """

    _lh: logging.Logger
    """
    Class logger handler
    """

    _cmd: Optional[List[str]]
    """Assembled Commandline"""

    _input_fasta: str
    """
    Input reference FASTA,
    can be DNA or transcript cDNA,
    can contain more than 1 entries.
    """

    _output_fastq_prefix: str
    """
    Prefix for output FASTQs.
    
    For LLRG that generates pair end (PE) reads,
    the generated read will named {_output_fastq_prefix}_1.fq and {_output_fastq_prefix}_2.fq.
    
    For LLRG that generates single-end (SE) reads,
    the generated read will named {_output_fastq_prefix}.fq
    """

    _depth: Union[int, float]
    """
    Sequencing depth or coverage.
    
    .. warning:: This is NOT final read count!
    """

    _tmp_dir: str
    """
    Simulator-based temp directory name.
    """

    # Following fields are left for LLRGs.
    _llrg_name: str
    """Class attribute, indicating name"""

    _require_integer_depth: bool
    """Class attribute, indicating whether the depts should be converted to ineger"""

    _capture_stdout: bool
    """Whether this simulator pours data into stdout"""

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: Union[int, float]
    ):
        super(BaseLLRGAdapter, self).__init__()
        if not hasattr(self, "_llrg_name"):
            raise TypeError
        if not hasattr(self, "_require_integer_depth"):
            raise TypeError
        if not hasattr(self, "_capture_stdout"):
            raise TypeError
        self._input_fasta = input_fasta
        self._output_fastq_prefix = output_fastq_prefix
        self._depth = int(depth) if self._require_integer_depth else depth
        self._lh = get_logger(__name__)
        self._cmd = None

    @abstractmethod
    def _pre_execution_hook(self) -> None:
        """
        Additional steps to do before starting the simulator process

        :raises LLRGException: If error occurs in LLRG-level.
        """
        raise NotImplementedError

    @abstractmethod
    def _rename_file_after_finish_hook(self):
        """
        Move the file into desired destination.
        May involve (de-)compression.

        :raises LLRGException: If error occurs in LLRG-level.
        """
        raise NotImplementedError

    def run(self):
        if self._cmd is None:
            self._lh.error("Commandline Assembly Failed!")
            return
        if not file_system.file_exists(self._input_fasta):
            self._lh.error("Commandline Assembly Failed!")
            return
        try:
            os.makedirs(self._tmp_dir, exist_ok=True)
        except (OSError, PermissionError, FileNotFoundError):
            self._lh.error(f"Failed to create temporary directory at %s", self._tmp_dir)
            return
        try:
            self._pre_execution_hook()
        except LLRGException as e:
            self._lh.error("Exception %s caught at pre-execution time", str(e))
            return

        subprocess_log_file_path = os.path.join(self._tmp_dir, "llrg.log")
        if self._capture_stdout:
            with open(subprocess_log_file_path, "wb") as subprocess_log_handler, \
                    open(self._output_fastq_prefix + ".fq", "wb") as stdout_handler:
                retv = self._exec_subprocess(
                    self._cmd,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_handler,
                    stderr=subprocess_log_handler
                )
        else:
            with open(subprocess_log_file_path, "wb") as subprocess_log_handler:
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

    def __repr__(self) -> str:
        return f"LLRG {self._llrg_name}: {self._input_fasta} -> {self._output_fastq_prefix} [{self._depth}]"

    def __str__(self):
        return repr(self)

    @property
    def llrg_name(self) -> str:
        return self._llrg_name

    @property
    def input_fasta(self) -> str:
        return self._input_fasta

    @property
    def output_fastq_prefix(self) -> str:
        return self._output_fastq_prefix

    @property
    def depth(self) -> Union[int, float]:
        return self._depth

    @abstractmethod
    @property
    def is_pair_end(self) -> bool:
        raise NotImplementedError
