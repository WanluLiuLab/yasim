# TODO
# New LLRG adapters should NOT simulate -> merge
# Should merge while simulating
# Can be implemented using callback functions

import subprocess
import threading
from abc import abstractmethod
from typing import Union, List, Callable

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger


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

    _before_cmd_hooks: List[Callable[[], None]]
    _after_cmd_hooks: List[Callable[[], None]]
    _on_success_hooks: List[Callable[[], None]]
    _on_failure_hooks: List[Callable[[], None]]

    # Following fields are left for LLRGs.
    _llrg_name: str
    _require_integer_depth: bool

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
        self.input_fasta = input_fasta
        self.output_fastq_prefix = output_fastq_prefix
        self.depth = int(depth) if self._require_integer_depth else depth
        self.exename = exename
        self.lh = get_logger(__name__)
        self._before_cmd_hooks = []
        self._after_cmd_hooks = []
        self._on_failure_hooks = []
        self._on_success_hooks = []
        self.other_args = other_args

    @abstractmethod
    def _assemble_cmd_hook(self) -> List[str]:
        """
        Assemble the command line into a list of string.
        """
        raise NotImplementedError

    @abstractmethod
    def _rename_file_after_finish_hook(self):
        """
        Move the file into desired destination.
        May involve (de-)compression.
        """
        raise NotImplementedError

    @abstractmethod
    def run(self):
        """
        The Default run method.
        Designed for LLRGs that does not pour generated sequences into STDOUT.
        """
        raise NotImplementedError

    def run_simulator_as_process_with_stdout_capturing(
            self
    ):
        """
        Use a subprocess to run LLRG. It performs following tasks:

        - Assembles CMD.
        - Execute the process.

        This method would redirect stdout into desired output location,
        so can only be used in single-end conditions.

        :return: The return value of top-level LLRG process. 0 for normal.
        """
        cmd = self._assemble_cmd_hook()
        self.lh.debug(f"Subprocess {' '.join(cmd)} START")
        with open(self.output_fastq_prefix + ".log", "wb") as subprocess_log_handler, \
                open(self.output_fastq_prefix + ".fq", "wb") as stdout_handler:
            p = subprocess.Popen(
                cmd,
                stdin=subprocess.DEVNULL,
                stdout=stdout_handler,
                stderr=subprocess_log_handler
            )
            retv = p.wait()
        if retv == 0:
            self.lh.debug(f"Subprocess {' '.join(cmd)} FIN")
            self._rename_file_after_finish_hook()
        else:
            self.lh.debug(f"Subprocess {' '.join(cmd)} ERR={retv}")
        return retv

    def run_simulator_as_process(self) -> int:
        """
        Use a subprocess to run LLRG.

        :return: The return value of top-level LLRG process. 0 for normal.
        """
        cmd = self._assemble_cmd_hook()
        self.lh.debug(f"Subprocess {' '.join(cmd)} START")
        with open(self.output_fastq_prefix + ".log", "wb") as subprocess_log_handler:
            p = subprocess.Popen(
                cmd,
                stdin=subprocess.DEVNULL,
                stdout=subprocess_log_handler,
                stderr=subprocess_log_handler
            )
            retv = p.wait()
        if retv == 0:
            self.lh.debug(f"Subprocess {' '.join(cmd)} FIN")
            self._rename_file_after_finish_hook()
        else:
            self.lh.debug(f"Subprocess {' '.join(cmd)} ERR={retv}")
        return retv
