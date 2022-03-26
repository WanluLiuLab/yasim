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
from typing import List, Union, Tuple, Dict, Optional, IO

from commonutils.io import determine_file_line_endings
from commonutils.io.file_system import file_exists
from commonutils.io.safe_io import get_reader, get_writer
from commonutils.io.tqdm_reader import get_tqdm_line_reader
from commonutils.stdlib_helper.logger_helper import chronolog, get_logger

lh = get_logger(__name__)

__all__ = [
    '_BaseFastaView',
    'FastaView'
]

QueryTuple = Union[Tuple[str, int, int], Tuple[str, int], Tuple[str]]


class _BaseFastaView:
    """
    Base class of other backends.
    """
    filename: str
    """
    Filename to read
    """

    full_header: bool
    """
    Whether to read in all header.
    If False (default), will read until seeing space or tab.
    See :py:mod:`pybedtools` for more details.
    """

    backend: str
    """
    The backend to use.
    """

    @chronolog(display_time=True)
    def __init__(self, filename: str, full_header: bool = False):
        self.full_header = full_header
        self.filename = filename
        self.backend = ''

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

    @abstractmethod
    def get_chr_length(self, chromosome: str) -> int:
        """
        Get length of specific chromosome.
        """
        pass

    @property
    @abstractmethod
    def chr_names(self) -> List[str]:
        """
        get chromosome names
        """
        pass

    def is_valid_region(self, chromosome: str, from_pos: int, to_pos: int):
        """
        Whether a region is valid. See :py:func:`sequence` for details.

        :raises ValueError: Raise this error if region is not valid.
        """
        if not chromosome in self.chr_names:
            raise ValueError(f"Chr {chromosome} not found")
        chr_len = self.get_chr_length(chromosome)
        if from_pos < 0 or from_pos > chr_len:
            raise ValueError(f"Seek {from_pos} too far")
        if to_pos != -1 and to_pos < 0 or to_pos > chr_len:
            raise ValueError(f"Seek {to_pos} too far")

    def __len__(self) -> int:
        return len(self.chr_names)

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
        """
        Write content of this FASTA view to file
        """
        with get_writer(output_filename) as writer:
            for k in self.chr_names:
                fa_str = f">{k}\n{self.sequence(k)}\n"
                writer.write(fa_str)

    def query(self, query: QueryTuple) -> str:
        """
        :py:func:`sequence` with another interface.
        """
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

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class _MemoryAccessFastaView(_BaseFastaView):
    """
    Fasta whose sequences are read into memory.
    Extremely fast but need lots of memory. Suitable for small files.
    """

    _all_dict: Dict[str, str]
    """
    Dict[chromosome_name, sequence]
    """

    @property
    def chr_names(self) -> List[str]:
        return list(self._all_dict.keys())

    def get_chr_length(self, chromosome: str) -> int:
        return len(self._all_dict[chromosome])

    def __init__(self, filename: str, full_header: bool = False):
        super().__init__(filename, full_header)
        self.backend = 'tetgs_in_memory'
        self._all_dict = {}  # For in-memory reader, will read in all sequences
        self._read_into_mem()

    @chronolog(display_time=True)
    def _read_into_mem(self) -> None:
        """
        Read FASTA into memory
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
                    seq = ''
                    line_len = 0
                if self.full_header:
                    chr_name = line[1:].strip()
                else:
                    chr_name = line[1:].strip().split(' ')[0].split('\t')[0]
            else:
                seq = seq + line
                if line_len == 0:
                    line_len = len(seq)
        if chr_name != '':
            self._all_dict[chr_name] = seq

    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1):
        self.is_valid_region(chromosome, from_pos, to_pos)
        if to_pos == -1:
            to_pos = self.get_chr_length(chromosome)
        return self._all_dict[chromosome][from_pos:to_pos]


class FastaIndexEntry:
    """
    An entry from ``.fai`` files
    """

    name: str
    """
    Chromosome Name
    """

    length: int
    """
    Chromosome length
    """

    offset: int
    """
    How far to `seek` to find this chromosome
    """

    line_blen: int
    """
    Line length in bytes without newline
    """

    line_len: int
    """
    Line length with newlines
    """

    @classmethod
    def from_fai_str(cls, fai_str: str):
        new_instance = cls()
        fields = fai_str.rstrip().split("\t")
        if len(fields) != 5:
            raise ValueError(f"Illegal record: {fai_str}. Need to have 5 fields.")
        new_instance.name = fields[0]
        new_instance.length = int(fields[1])
        new_instance.offset = int(fields[2])
        new_instance.line_blen = int(fields[3])
        new_instance.line_len = int(fields[4])
        return new_instance

    def __repr__(self):
        return "\t".join((
            self.name,
            str(self.length),
            str(self.offset),
            str(self.line_blen),
            str(self.line_len)
        ))

    def __str__(self):
        return repr(self)


FAI_INDEX_TYPE = Dict[str, FastaIndexEntry]


def create_fai(
        filename: str,
        index_filename: str
) -> FAI_INDEX_TYPE:
    """
    Create an FAI index for Fasta

    Do not use this feature on full headers.
    """
    return_fai: FAI_INDEX_TYPE = {}
    name = ""
    offset = 0
    line_len = 0
    line_blen = 0
    length = 0

    def append():
        new_fai_index = FastaIndexEntry()
        new_fai_index.name = name
        new_fai_index.line_len = line_len
        new_fai_index.line_blen = line_blen
        new_fai_index.offset = offset
        new_fai_index.length = length
        return_fai[name] = new_fai_index

    with get_tqdm_line_reader(filename, newline=determine_file_line_endings(filename)) as reader:
        while True:
            line = reader.readline()
            if not line:
                break
            if line == '':
                continue
            if line[0] == '>':  # FASTA header
                if name != '':
                    append()
                    offset = 0
                    line_len = 0
                    line_blen = 0
                    length = 0
                name = line[1:].strip().split(' ')[0].split('\t')[0]
                offset = reader.tell()
            else:
                if line_len == 0:
                    line_len = len(line)
                if line_blen == 0:
                    line_blen = len(line.rstrip('\r\n'))
                length += len(line.rstrip('\r\n'))
        if name != '':
            append()
    with get_writer(index_filename) as writer:
        for fai_record in return_fai.values():
            writer.write(str(fai_record) + "\n")
    return return_fai


class _DiskAccessFastaView(_BaseFastaView):
    """
    Fasta whose sequence is NOT read into memory, with :py:mod:``tetgs`` backend.
    Slow but memory-efficient.

    # FIXME: Error handling one-line FASTA
    """

    _fd: IO
    """
    Underlying file descriptor
    """

    _fai: FAI_INDEX_TYPE

    def get_chr_length(self, chromosome: str) -> int:
        return self._fai[chromosome].length

    @property
    def chr_names(self) -> List[str]:
        return list(self._fai.keys())

    def __init__(self, filename: str, full_header: bool = False):
        super().__init__(filename, full_header)
        self.backend = 'tetgs'
        # If has prebuilt index file, read it
        self._fd = get_reader(self.filename)
        index_filename = self.filename + ".fai"
        if not file_exists(index_filename) or \
                os.path.getmtime(index_filename) - os.path.getmtime(filename) < 0:
            create_fai(self.filename, index_filename)
        self._read_index_from_fai(index_filename)

    @chronolog(display_time=True)
    def _read_index_from_fai(self, index_filename: str):
        self._fai = {}
        for line in get_tqdm_line_reader(index_filename):
            index_entry = FastaIndexEntry.from_fai_str(line)
            self._fai[index_entry.name] = index_entry

    def sequence(self, chromosome: str, from_pos: int = 0, to_pos: int = -1) -> str:
        self.is_valid_region(chromosome, from_pos, to_pos)
        chr_fai = self._fai[chromosome]
        """FAI record of this chromosome"""

        if to_pos == -1:
            to_pos = self.get_chr_length(chromosome)

        len_newline = chr_fai.line_len - chr_fai.line_blen
        """Length of newlines"""

        self._fd.seek(chr_fai.offset + from_pos // chr_fai.line_blen * len_newline + from_pos)
        prev_resid = from_pos % chr_fai.line_blen
        """Previous residue. Where the reader is on"""

        lines_to_read = (to_pos - from_pos + prev_resid) // chr_fai.line_blen
        """
        How many full-length line to read
        """

        if (to_pos % chr_fai.line_blen) == 0:
            lines_to_read -= 1
        rets = self._fd.read(
            to_pos - from_pos + lines_to_read
        ).replace('\n', '').replace('\r', '')
        return rets

    def close(self):
        try:
            self._fd.close()
        except AttributeError:
            pass


class FastaView(_BaseFastaView, ABC):
    """
    The major Fasta handler class, supporting multiple backends.
    """

    def __new__(cls,
                filename: str,
                full_header: bool = False,
                read_into_memory: Optional[bool] = None):
        """
        Initialize a _DiskFasta interface using multiple backends.

        :param filename: The file you wish to open.
        :param full_header: Whether to read full headers.
        :param read_into_memory: Whether to read into memory.
        """
        # return _MemoryAccessFastaView(filename, full_header)
        if read_into_memory:
            return _MemoryAccessFastaView(filename, full_header)
        else:
            return _DiskAccessFastaView(filename, full_header)
