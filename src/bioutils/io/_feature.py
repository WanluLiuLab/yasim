from typing import Iterator

from bioutils.datastructure import Feature
from bioutils.datastructure import Gff3Record
from bioutils.datastructure import GtfRecord
from commonutils import ioctl
from commonutils.tqdm_utils import tqdm_line_reader


class _FeatureIterator:
    filename:str = ""
    filetype:str = "UNKNOWN"
    def __init__(self, filename:str):
        self.filename = ioctl.ensure_input_existence(filename)

    def __iter__(self) -> Iterator[Feature]:
        pass

    def __repr__(self):
        return f"{self.filetype} Iterator for {self.filename}"

    __str__ = __repr__


class GtfIterator(_FeatureIterator):
    filetype: str = "GTF"
    def __iter__(self) -> Iterator[GtfRecord]:
        with tqdm_line_reader(self.filename) as reader:
            while True:
                line = reader.readline()
                if not line:
                    break
                line = line.strip()
                if line.startswith('#') or line == '':
                    continue
                yield GtfRecord.from_string(line)

class Gff3Iterator(_FeatureIterator):
    filetype: str = "GFF3"
    def __iter__(self) -> Iterator[Gff3Record]:
        with tqdm_line_reader(self.filename) as reader:
            while True:
                line = reader.readline()
                if not line:
                    break
                line = line.strip()
                if line.startswith('#') or line == '':
                    continue
                yield Gff3Record.from_string(line)


