"""
fasta_view.py -- General FASTA reader

Can provide random access to FASTA files, compressed or non-compressed.

Highlights: This utility can read all format supported by :py:mod:`commonutils.io`, while others require Block GZipped ones.

.. note::
    Although this module supports all format supported by :py:mod:`commonutils.io`,
    it is recommended for user to compress their files using ``bgzip`` and index them using ``tabix``.

.. warning::
    This module uses 0-based ``[)`` indexing!
"""

import os
from abc import abstractmethod, ABC
from typing import List, Union, Tuple, Dict

from commonutils.io.safe_io import get_reader, get_writer
from commonutils.io.tqdm_reader import get_tqdm_line_reader
from commonutils.stdlib_helper.logger_helper import chronolog, get_logger

lh = get_logger(__name__)

__all__ = [
    '_FastaView',
    'FastaView',
    '__version__'
]

__version__ = 0.2

QueryTuple = Union[Tuple[str, int, int], Tuple[str, int], Tuple[str]]


class _FastaView:
    """
    Base class of other backends.
    """
    filename: str
    """
    Filename to read
    """

    relname: str
    """
    The absolute path of the filename.
    """

    all_header: bool
    """
    Whether to read in all header.
    If False (default), will read until seeing space or tab.
    See :py:mod:`pybedtools` for more details.
    """

    backend: str
    """
    The backend to use.
    """

    chr_dict: Dict[str, Tuple[int, int, int]]
    """ 
    Format: ``{chr: (length, seek_start, line_width)}``
    The line_width is the 4th (last but one) field of .fai files
    """

    @chronolog(display_time=True)
    def __init__(self, filename: str, all_header: bool = False):
        self.all_header = all_header
        self.filename = filename
        self.backend = ''
        self._chr_dict = dict()

    @abstractmethod
    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> str:
        """
        Get sequence from FASTA with 0-based [) indexes

        To read until end, use -1.

        :param chromosome: Chromosome name
        :param from_pos: From which position
        :param to_pos: To which position, use -1 for end
        """
        pass

    def is_valid_region(self, chromosome: str, from_pos: int, to_pos: int):
        """
        Whether a region is valid. See :py:func:`sequence` for details.

        :raises ValueError: Raise this error if region is not valid.
        """
        try:
            chr_len = self._chr_dict[chromosome][0]
        except KeyError as e:
            raise ValueError(f"Chr {chromosome} not found") from e
        if from_pos < 0 or from_pos > chr_len:
            raise ValueError(f"Seek {from_pos} too far")
        if to_pos != -1 and to_pos < 0 or to_pos > chr_len:
            raise ValueError(f"Seek {to_pos} too far")

    def __len__(self) -> int:
        return len(self._chr_dict)

    def __repr__(self):
        try:
            return f"Fasta from {self.filename} with backend {self.backend}, len={len(self)}"
        except AttributeError:
            return "Fasta being constructed"

    __str__ = __repr__

    def close(self):
        """
        Safely close a Fasta.
        """
        pass

    def __del__(self):
        self.close()

    def to_file(self, output_filename: str):
        with get_writer(output_filename) as writer:
            for k in self._chr_dict.keys():
                fa_str = f">{k}\n{self.sequence(k)}\n"
                writer.write(fa_str)

    def query(self, query: QueryTuple) -> str:
        return self.sequence(*query)

    def subset(self, output_filename: str, querys: List[QueryTuple]):
        with get_writer(output_filename) as writer:
            for query in querys:
                fa_str = f">{query[0]}\n{self.query(query)}\n"
                writer.write(fa_str)

    def subset_chr(self, output_filename: str, querys: List[str]):
        with get_writer(output_filename) as writer:
            for query in querys:
                fa_str = f">{query}\n{self.sequence(query)}\n"
                writer.write(fa_str)


class _MemoryAccessFastaView(_FastaView):
    """
    Fasta whose sequences are read into memory.
    Extremely fast but need lots of memory. Suitable for small files.
    """

    def __init__(self, filename: str, all_header: bool = False):
        super().__init__(filename, all_header)
        self.backend = 'tetgs_in_memory'
        self._all_dict = dict()  # For in-memory reader, will read in all sequences
        self._read_into_mem()

    @chronolog(display_time=True)
    def _read_into_mem(self) -> None:
        """
        Construct _chr_dict from file and read into memory
        """
        chr_name = ""
        seq = ""
        line_len = 0
        for line in get_tqdm_line_reader(self.filename):
            if line == "":
                continue
            if line[0] == '>':  # FASTA header
                if chr_name != '':
                    self._all_dict[chr_name] = seq
                    self._chr_dict[chr_name][0] = len(seq)
                    self._chr_dict[chr_name][2] = line_len
                    seq = ''
                    line_len = 0
                if self.all_header:
                    chr_name = line[1:].strip()
                else:
                    chr_name = line[1:].strip().split(' ')[0].split('\t')[0]
                self._chr_dict[chr_name] = [0, 0, 0]
            else:
                seq = seq + line
                if line_len == 0:
                    line_len = len(seq)
        if chr_name != '':
            self._all_dict[chr_name] = seq
            self._chr_dict[chr_name][0] = len(seq)
            self._chr_dict[chr_name][2] = line_len

    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1):
        self.is_valid_region(chromosome, from_pos, to_pos)
        if to_pos == -1:
            to_pos = self._chr_dict[chromosome][0]
        return self._all_dict[chromosome][from_pos:to_pos]


class _DiskAccessFastaView(_FastaView):
    """
    Fasta whose sequence is NOT read into memory, with :py:mod:``tetgs`` backend.
    Slow but memory-efficient.

    # FIXME: Error handling one-line FASTA
    """

    def __init__(self, filename: str, all_header: bool = False):
        super().__init__(filename, all_header)
        self.backend = 'tetgs'
        # If has prebuilt index file, read it
        self._reader = get_reader(self.filename)
        if os.path.exists(self.filename + '.fai'):
            self._read_index_from_fai(self.filename + '.fai')
        else:
            self._create_index()

    @chronolog(display_time=True)
    def _read_index_from_fai(self, index_filename: str):
        with get_reader(index_filename) as fh:
            while True:
                line = fh.readline()
                if not line:
                    break
                lines = line.rstrip().split('\t')
                self._chr_dict[lines[0]] = [int(lines[1]), int(lines[2]), int(lines[3])]

    @chronolog(display_time=True)
    def _create_index(self):
        """
        Construct _chr_dict from file without read into memory
        """
        chr_name = ""
        seq_len = 0
        line_len = 0
        with get_tqdm_line_reader(self.filename) as reader:  # do not use get_tqdm_line_reader; it does not `tell`.
            while True:
                line = reader.readline()
                if not line:
                    break
                line = line.rstrip()
                if line == '':
                    continue
                if line[0] == '>':  # FASTA header
                    if chr_name != '':
                        self._chr_dict[chr_name][0] = seq_len
                        if line_len == 0:
                            line_len = seq_len
                        self._chr_dict[chr_name][2] = line_len
                        seq_len = 0
                        line_len = 0
                    if self.all_header:
                        chr_name = line[1:].strip()
                    else:
                        chr_name = line[1:].strip().split(' ')[0].split('\t')[0]
                    self._chr_dict[chr_name] = [0, reader.tell(), 0]
                else:
                    if line_len == 0:
                        line_len = seq_len
                    seq_len = seq_len + len(line)
            if chr_name != '':
                self._chr_dict[chr_name][0] = seq_len
                if line_len == 0:
                    line_len = seq_len
                self._chr_dict[chr_name][2] = line_len

    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> str:
        self.is_valid_region(chromosome, from_pos, to_pos)
        if to_pos == -1:
            to_pos = self._chr_dict[chromosome][0]
        line_len = self._chr_dict[chromosome][2]
        self._reader.seek(self._chr_dict[chromosome][1] + from_pos // line_len + from_pos)
        prev_resid = from_pos % line_len
        """
        Previous residue. Where the reader is on.
        """

        lines_to_read = (to_pos - from_pos + prev_resid) // line_len
        """
        How many full-length line to read
        """

        if (to_pos % line_len) == 0:
            lines_to_read -= 1
        rets = self._reader.read(to_pos - from_pos + lines_to_read).replace('\n', '')
        return rets

    def close(self):
        try:
            self._reader.close()
        except AttributeError:
            pass


class FastaView(_FastaView, ABC):
    """
    The major Fasta handler class, supporting multiple backends.
    """

    def __new__(cls,
                filename: str,
                all_header: bool = False,
                read_into_memory: bool = False):
        """
        Initialize a _DiskFasta interface using multiple backends.

        :param filename: The file you wish to open.
        :param all_header: Whether to read full headers.
        :param read_into_memory: Whether to read into memory.
        """
        return _MemoryAccessFastaView(filename, all_header)
        # if read_into_memory:
        #     return _MemoryAccessFastaView(filename, all_header)
        # else:
        #     return _DiskAccessFastaView(filename, all_header)
