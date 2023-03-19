"""
yasim.llrg_adapter -- Adapters for LLRGs

Here contains LLRG adapters that should have been used in DNA-Seq.
"""

from __future__ import annotations

import logging
import os.path
import subprocess
import threading
from abc import abstractmethod
from typing import Union, List, IO, Optional, Iterable

from labw_utils.commonutils.io.safe_io import get_reader, get_writer, file_system
from labw_utils.commonutils.shell_utils import wc_c
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.llrg import enhanced_which

COPY_BUFSIZE = 1024 * 1024 if os.name == 'nt' else 64 * 1024
"""Copy buffer size from shutil."""


class LLRGException(RuntimeError):
    """Some error raised due to (improperly configured) LLRG"""

    _contents: str

    def __init__(self, contents: str):
        """
        Default Initializer.
        :param contents: What to be put in logs.
        """
        self._contents = contents

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return self._contents


class NoOutputFileException(LLRGException):
    """Error raised if no output file was found even when LLRG exited normally."""
    ...


class EmptyOutputFileException(LLRGException):
    """Error raised if empty output file was found even when LLRG exited normally."""
    ...


class LLRGFailException(LLRGException):
    """Error raised LLRG exited abnormally."""
    ...


class LLRGInitializationException(LLRGException):
    """Error raised in initialization or pre-execution stage."""
    ...


def enhanced_copyfileobj(src_fd: IO, dst_fd: IO) -> int:
    """
    :py:func:`shutil.copyfileobj` with number of bytes copied.
    """
    c_len = 0
    while True:
        buf = src_fd.read(COPY_BUFSIZE)
        if not buf:
            break
        dst_fd.write(buf)
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
    - Move generated files to ``_dst_fastq_file_prefix``.
    """

    _src_fasta_file_path: str
    _dst_fastq_file_prefix: str
    _depth: Union[int, float]

    # Following fields are left for LLRGs.

    _lh: logging.Logger
    """
    Class logger handler
    """

    _cmd: Optional[List[str]]
    """Assembled Commandline"""

    _tmp_dir: str
    """
    Simulator-based temp directory name.
    """

    _exception: Optional[LLRGException]
    """Exception raised during initialization, pre-execution, execution or post-execution of LLRG."""

    # Following fields are left for LLRG class variables..
    llrg_name: str
    """Class attribute, indicating name. Should be Final."""

    _require_integer_depth: bool
    """Class attribute, indicating whether the depth should be converted to integer. Should be Final."""

    _capture_stdout: bool
    """Whether this simulator pours data into stdout. Should be Final."""

    def __init__(
            self,
            src_fasta_file_path: str,
            dst_fastq_file_prefix: str,
            depth: Union[int, float],
            llrg_executable_path: str,
            is_trusted: bool,
            **kwargs
    ):
        """
        Initializer.

        :param src_fasta_file_path: Path of source FASTA.
        :param dst_fastq_file_prefix: Prefix od destination FASTQ.
        :param depth: Targeted sequencing depth. Would NOT be related to actual sequencing depth!
        :param llrg_executable_path: Path to LLRG Executable.
        :param kwargs: Other miscellaneous arguments.
        :param is_trusted: Whether to skip input validation test.
        :raise LLRGInitializationException: On error.
        """
        # To developers: This function should raise LLRGInitializationException only!
        super(BaseLLRGAdapter, self).__init__()

        # Check for whether class variables have been properly set.
        if not hasattr(self, "llrg_name"):
            raise TypeError
        if not hasattr(self, "_require_integer_depth"):
            raise TypeError
        if not hasattr(self, "_capture_stdout"):
            raise TypeError
        self._src_fasta_file_path = src_fasta_file_path
        self._dst_fastq_file_prefix = os.path.abspath(dst_fastq_file_prefix)
        self._depth = int(depth) if self._require_integer_depth else depth
        self._lh = get_logger(__name__)
        self._cmd = None
        self._tmp_dir = self._dst_fastq_file_prefix + ".tmp.d"
        self._exception = None

        # Validate input.
        if not is_trusted:
            if self._depth <= 0:
                raise LLRGInitializationException(f"Depth {self._depth} too low")
            if not file_system.file_exists(self._src_fasta_file_path):
                raise LLRGInitializationException(f"FASTA {self._src_fasta_file_path} not found!")
            try:
                os.makedirs(self._tmp_dir, exist_ok=True)
            except (OSError, PermissionError, FileNotFoundError) as e:
                raise LLRGInitializationException("MKTEMP Failed!") from e
            try:
                _ = enhanced_which(llrg_executable_path)
            except FileNotFoundError as e:
                raise LLRGInitializationException(
                    f"LLRG Executable {llrg_executable_path} not found or not executable!"
                ) from e

    @property
    def exception(self) -> str:
        """
        Exception raised during initialization, pre-execution, execution or post-execution of LLRG.

        :return: Following values:

            - ``NORMAL``: If no exception occurs
            - ``EmptyOutFile``: If LLRG exited normally but with empty output file.
            - ``NoOutputFile``: If LLRG exited normally but with no output file.
            - ``LLRGFail``: If LLRG exited abnormally.
            - ``InitFail``: if pre-execution of initialization hook failed.
            - ``UNKNOWN``: Other errors.
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
        Additional steps to take before starting the simulator process

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

    def _run_llrg_hook(self) -> None:
        """
        Execute LLRG as a process

        :raise LLRGFailException: On running failures.
        """
        subprocess_log_file_path = os.path.join(self._tmp_dir, "llrg.log")
        dst_fastq_file_path = self._dst_fastq_file_prefix + ".fq"
        if self._capture_stdout:
            with get_writer(subprocess_log_file_path, is_binary=True) as subprocess_log_handler, \
                    get_writer(dst_fastq_file_path, is_binary=True) as stdout_handler:
                retv = self._exec_subprocess(
                    self._cmd,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_handler,
                    stderr=subprocess_log_handler
                )
            if wc_c(dst_fastq_file_path) < 2:
                raise EmptyOutputFileException(f"Output {dst_fastq_file_path} empty!")
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

    def run(self) -> None:
        """Execute LLRG"""
        try:
            if self._cmd is None:
                raise LLRGInitializationException("Commandline Assembly Failed!")
            self._pre_execution_hook()
        except LLRGException as e:
            self._lh.warning("Exception %s caught at pre-execution time", str(e))
            self._exception = e
            return
        try:
            self._run_llrg_hook()
        except LLRGException as e:
            self._lh.warning("Exception %s caught at execution time", str(e))
            self._exception = e
            return
        try:
            self._post_execution_hook()
        except LLRGException as e:
            self._lh.warning("Exception %s caught at post-execution hook", str(e))
            self._exception = e
            return

    def _exec_subprocess(
            self,
            cmd: List[str],
            stdin: Union[IO, int],
            stdout: Union[IO, int],
            stderr: Union[IO, int]
    ) -> int:
        """
        Wrapper of :py:class:`subprocess.Popen` with logs.
        Will not raise exception if error occurs.

        :param cmd: Commandline arguments.
        :param stdin: File descriptor for STDIN.
        :param stdout: File descriptor for STDOUT.
        :param stderr: File descriptor for STDERR.
        :return: Return value of the process.
        """
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
            self._lh.warning(f"Subprocess {' '.join(cmd)} ERR={retv}")
        return retv

    def __repr__(self) -> str:
        return f"LLRG {self.llrg_name}: {self._src_fasta_file_path} -> {self._dst_fastq_file_prefix} [{self._depth}]"

    def __str__(self):
        return repr(self)

    @property
    def input_fasta(self) -> str:
        """
        Input reference FASTA.

        - DNA: Can contain more than one entry.
        - Isoform cDNA: Should contain only one entry.
        """
        return self._src_fasta_file_path

    @property
    def output_fastq_prefix(self) -> str:
        """
        Prefix for output FASTQs.

        For LLRG that generates pair end (PE) reads,
        the generated read will be named ``{_dst_fastq_file_prefix}_1.fq`` and ``{_dst_fastq_file_prefix}_2.fq``.

        For LLRG that generates single-end (SE) reads,
        the generated read will be named ``{_dst_fastq_file_prefix}.fq``
        """
        return self._dst_fastq_file_prefix

    @property
    def depth(self) -> Union[int, float]:
        """
        Sequencing depth or coverage.

        .. warning:: This is NOT final read count!
        """
        return self._depth

    @property
    @abstractmethod
    def is_pair_end(self) -> bool:
        """
        Whether this simulator is pair-end.
        """
        raise NotImplementedError
