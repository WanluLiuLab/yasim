# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: new_gtf.py -- General-purposed GTF reader.
#
#  VERSION HISTORY:
#  VERSION HISTORY:
#  2021-08-16 0.1  : Purposed and added by YU Zhejian.
#  2021-08-24 0.1  : IntervalTree added.
#  2021-08-24 0.2  : Refactored as is designed by YUAN Ruihong.
#  2021-08-30 0.3  : Refactored as is designed by YUAN Ruihong.
#  2021-09-12 0.4  : Cythonized.
#
# ==============================================================================

"""
gtf.py -- General-purposed GTF reader.

Just as Fasta, this general-purposed reader should support multiple backends.
But currently it supports tetgs only.

.. warning::
    This file uses 1-based [] indexing!
"""
import time
from collections import OrderedDict
from typing import Iterator, Union, Optional, Iterable, List

import intervaltree

from bioutils.io import fasta
from commonutils import ioctl, logger
from commonutils.tqdm_importer import tqdm

try:
    from bioutils.io.gtf._c_gtf_record import GtfRecord
except ImportError:
    from bioutils.io.gtf._py_gtf_record import GtfRecord

from commonutils.tqdm_utils import tqdm_line_reader

lh = logger.get_logger(__name__)

__version__ = 0.4


class _GtfView:

    def __init__(self):
        self.attached_fasta = None

    def attach_external_fasta(self, fasta_filename: str):
        self.attached_fasta = fasta.FastaView(fasta_filename)

    def close_external_fasta(self):
        try:
            self.attached_fasta.close()
            self.attached_fasta = None
        except AttributeError:
            self.attached_fasta = None

    def __iter__(self) -> Iterator[GtfRecord]:
        pass

    def __del__(self):
        self.close_external_fasta()


class GtfIterator(_GtfView):

    def __init__(self, gtf_filename: str):
        super().__init__()
        self.gtf_filename = gtf_filename

    def __iter__(self) -> Iterator[GtfRecord]:
        realname = ioctl.ensure_input_existence(self.gtf_filename)

        with tqdm_line_reader(realname) as reader:
            while True:
                line = reader.readline()
                if not line:
                    break
                line = line.rstrip()
                if line.startswith('#') or line == '':
                    continue
                gtf_record = GtfRecord.from_string(line)  # 91.7% of time
                gtf_record.gtf_handler = self
                yield gtf_record  # TODO: Extremely slow for memory fragmentation


class SimpleGtfView(_GtfView):
    _comments: List[str]

    def __init__(self, gtf_filename_or_iterable: Optional[Union[str, Iterable[GtfRecord]]] = None):
        super().__init__()
        if gtf_filename_or_iterable and isinstance(gtf_filename_or_iterable, str):
            realname = ioctl.ensure_input_existence(gtf_filename_or_iterable)
            self.__init__(GtfIterator(realname))
            return

        self._comments = []  # of strings
        self._contigs = OrderedDict()  # of interval trees

        if gtf_filename_or_iterable and isinstance(gtf_filename_or_iterable, Iterable):
            for gtf_record in gtf_filename_or_iterable:  # 45.4% of time
                gtf_record.gtf_handler = self
                self.add_record(gtf_record)  # 54.6% of time

    @property
    def contig_names(self):
        return self._contigs.keys()

    # @profile
    def fetch(self,
              seqname: str,
              start: int = 0,
              end: Optional[int] = None
              ) -> Iterable[GtfRecord]:
        """
        Get an iterator of :py:class:`GtfRecord` overlapping with given interval.

        :param seqname: The seqname to search.
        :param start: 1-based inclusive interval start.
        :param end: 1-based inclusive interval end.
        :return: Iterator of overlapping record.
        """
        if seqname not in self._contigs.keys():
            return
        if end == None:
            interval_tree = self._contigs[seqname]
        else:
            interval_tree = self._contigs[seqname][start:end + 1]
        for interval in interval_tree:
            yield interval.data

    def add_record(self, record: GtfRecord):
        """
        Add a new GTF record

        :param record: New record to be added.
        """
        record.gtf_handler = self
        if record.seqname not in self._contigs.keys():
            self._contigs[record.seqname] = intervaltree.IntervalTree()
        self._contigs[record.seqname].addi(record.start, record.end, record)  # 98.9% of time

    def to_file(self, filename: str):
        """
        Write the contents to another file.

        :param filename: Filename to write.
        """
        self._comments.append(f'# Time {time.asctime()}\n')
        realname = ioctl.ensure_output_existence(filename)
        with ioctl.get_writer(realname) as writer:
            for item in tqdm(iterable=self._comments, desc=f"Writing to {filename}: Step 1 of 2 - comments..."):
                writer.write(str(item) + '\n')
            for item in tqdm(iterable=self, desc=f"Writing to {filename}: Step 2 of 2 - records...",
                             total=len(self)):
                writer.write(str(item) + '\n')

    def __iter__(self) -> Iterator[GtfRecord]:
        """
        :return: An iterator of all records in order.
        """
        for contig in self._contigs.values():
            for interval in contig:
                yield interval.data

    def __len__(self) -> int:
        """
        Get number of records.

        :return: Number of records.
        """
        return sum(len(itree) for itree in self._contigs.values())
