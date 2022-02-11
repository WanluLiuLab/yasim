from typing import Iterator, Union, Optional, List

from bioutils.datastructure import Feature
from bioutils.datastructure.gff_gtf_record import Gff3Record
from bioutils.datastructure.gff_gtf_record import GtfRecord
from commonutils import ioctl
from commonutils.tqdm_utils import tqdm_line_reader


class _FeatureIterator:
    filename: str = ""
    filetype: str = "UNKNOWN"

    def __init__(self, filename: str):
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


class _FeatureWriter:
    @staticmethod
    def write(
            iterable: Union[_FeatureIterator, Iterator[Feature]],
            output_filename: str,
            prefix_annotations: Optional[List[str]] = None
    ):
        with ioctl.get_writer(output_filename) as writer:
            if prefix_annotations is not None:
                for annotation in prefix_annotations:
                    writer.write('#' + annotation + '\n')
            for feature in iterable:
                writer.write(str(feature) + '\n')


class GtfWriter(_FeatureWriter):
    pass


class Gff3Writer(_FeatureWriter):
    pass
