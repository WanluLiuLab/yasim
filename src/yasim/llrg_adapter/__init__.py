import logging
import os.path
import subprocess
import threading
from abc import abstractmethod
from typing import Union, List, IO, Optional, Iterable

from labw_utils.commonutils.io import get_reader, get_writer, file_system
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

COPY_BUFSIZE = 1024 * 1024 if os.name == 'nt' else 64 * 1024
"""Copy buffer size from shutil."""


class LLRGException(RuntimeError):
    """Some error raised due to (improperly configured) LLRG"""

    def __init__(self, contents: str):
        self.contents = contents

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return self.contents


class NoOutputFileException(LLRGException):
    ...


class EmptyOutputFileException(LLRGException):
    ...


class LLRGFailException(LLRGException):
    ...


class LLRGInitializationException(LLRGException):
    """Error raised in initialization or pre-execution stage."""
    ...


def enhanced_copyfileobj(fsrc: IO, fdst: IO, length: int = 0) -> int:
    """
    :py:func:`shutil.copyfileobj` with number of bytes copied.
    """
    c_len = 0
    if not length:
        length = COPY_BUFSIZE
    fsrc_read = fsrc.read
    fdst_write = fdst.write
    while True:
        buf = fsrc_read(length)
        if not buf:
            break
        fdst_write(buf)
        c_len += len(buf)
    return c_len


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
            c_len = enhanced_copyfileobj(r1, w1)
    except FileNotFoundError as e:
        raise NoOutputFileException(f"Copy file {in_fn} not found!") from e
    except (OSError, PermissionError, IOError) as e:
        raise LLRGException(f"Copy file {in_fn} -> {out_fn} failed!") from e
    if c_len == 0:
        raise EmptyOutputFileException(f"Copy file {in_fn} empty!")


def automerge(in_fns: Iterable[str], out_fn: str) -> None:
    """
    Merge multiple files into one. See :py:func:`autocopy`.
    """
    c_len = 0
    with get_writer(out_fn, is_binary=True) as w1:
        for in_fn in in_fns:
            try:
                with get_reader(in_fn, is_binary=True) as r1:
                    c_len += enhanced_copyfileobj(r1, w1)
            except FileNotFoundError as e:
                raise NoOutputFileException(f"Copy file {in_fn} not found!") from e
            except (OSError, PermissionError, IOError) as e:
                raise LLRGException(f"Copy file {in_fn} -> {out_fn} failed!") from e
    if c_len == 0:
        raise EmptyOutputFileException(f"Copy file {in_fn} empty!")


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

    _exception: Optional[LLRGException]

    def __init__(
            self,
            input_fasta: str,
            output_fastq_prefix: str,
            depth: Union[int, float],
            **kwargs
    ):
        # To developers: This function should raise LLRGInitializationException only!
        super(BaseLLRGAdapter, self).__init__()
        if not hasattr(self, "_llrg_name"):
            raise TypeError
        if not hasattr(self, "_require_integer_depth"):
            raise TypeError
        if not hasattr(self, "_capture_stdout"):
            raise TypeError
        self._input_fasta = input_fasta
        self._output_fastq_prefix = os.path.abspath(output_fastq_prefix)
        self._depth = int(depth) if self._require_integer_depth else depth
        self._lh = get_logger(__name__)
        self._cmd = None
        self._tmp_dir = self._output_fastq_prefix + ".tmp.d"
        self._exception = None

    @property
    def exception(self) -> str:
        """
        Get Exception status.

        :return:
        """
        if self._exception is None:
            return "NORMAL"
        elif isinstance(self._exception, EmptyOutputFileException):
            return "EmptyOutFile"
        elif isinstance(self._exception, NoOutputFileException):
            return "NoOutputFile"
        elif isinstance(self._exception, LLRGFailException):
            return "LLRGFail"
        elif isinstance(self._exception, LLRGInitializationException):
            return "InitFail"
        else:
            return "UNKNOWN"

    @abstractmethod
    def _pre_execution_hook(self) -> None:
        """
        Additional steps to do before starting the simulator process

        :raises LLRGException: If error occurs in LLRG-level.
        """
        raise NotImplementedError

    @abstractmethod
    def _post_execution_hook(self):
        """
        Move the file into desired destination.
        May involve (de-)compression.

        :raises LLRGException: If error occurs in LLRG-level.
        """
        raise NotImplementedError

    def _llrg_initialization_hook(self) -> None:
        """
        Perform pre-flight check for LLRG.

        :raise LLRGInitializationException: On initiation failures.
        """
        if self._cmd is None:
            raise LLRGInitializationException("Commandline Assembly Failed!")
        if not file_system.file_exists(self._input_fasta):
            raise LLRGInitializationException(f"FASTA {self._input_fasta} not found!")
        try:
            os.makedirs(self._tmp_dir, exist_ok=True)
        except (OSError, PermissionError, FileNotFoundError) as e:
            raise LLRGInitializationException(f"MKTEMP Failed!") from e

    def _run_llrg_hook(self) -> None:
        """
        Execute LLRG as a process

        :raise LLRGFailException: On running failures.
        """
        subprocess_log_file_path = os.path.join(self._tmp_dir, "llrg.log")
        if self._capture_stdout:
            with get_writer(subprocess_log_file_path, is_binary=True) as subprocess_log_handler, \
                    get_writer(self._output_fastq_prefix + ".fq", is_binary=True) as stdout_handler:
                retv = self._exec_subprocess(
                    self._cmd,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_handler,
                    stderr=subprocess_log_handler
                )
        else:
            with get_writer(subprocess_log_file_path, is_binary=True) as subprocess_log_handler:
                retv = self._exec_subprocess(
                    self._cmd,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess_log_handler,
                    stderr=subprocess_log_handler
                )
        if retv != 0:
            raise LLRGFailException(f"Return value LLRG ({retv}) != 0")

    def run(self):
        try:
            self._llrg_initialization_hook()
            self._pre_execution_hook()
        except LLRGException as e:
            self._lh.error("Exception %s caught at pre-execution time", str(e))
            self._exception = e
            return
        try:
            self._run_llrg_hook()
        except LLRGException as e:
            self._lh.error("Exception %s caught at execution time", str(e))
            self._exception = e
            return
        try:
            self._post_execution_hook()
        except LLRGException as e:
            self._lh.error("Exception %s caught at post-execution hook", str(e))
            self._exception = e
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

    @property
    @abstractmethod
    def is_pair_end(self) -> bool:
        raise NotImplementedError
