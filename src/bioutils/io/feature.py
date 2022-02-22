from abc import abstractmethod
from typing import Iterator, Union, Optional, List, TextIO

from bioutils.datastructure import Feature
from bioutils.datastructure.gff_gtf_record import Gff3Record
from bioutils.datastructure.gff_gtf_record import GtfRecord
from commonutils import ioctl
from commonutils.tqdm_utils import tqdm_line_iterator


class _FeatureIterator:
    filename: str = ""
    filetype: str = "UNKNOWN"

    def __init__(self, filename: str):
        self.filename = ioctl.ensure_input_existence(filename)

    @abstractmethod
    def __iter__(self) -> Iterator[Feature]:
        pass

    def __repr__(self):
        return f"{self.filetype} Iterator for {self.filename}"

    __str__ = __repr__

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        """
        pass


class GtfIterator(_FeatureIterator):
    filetype: str = "GTF"

    def __iter__(self) -> Iterator[GtfRecord]:
        for line in tqdm_line_iterator(self.filename):
            if line.startswith('#') or line == '':
                continue
            yield GtfRecord.from_string(line)


class Gff3Iterator(_FeatureIterator):
    filetype: str = "GFF3"

    def __iter__(self) -> Iterator[Gff3Record]:
        for line in tqdm_line_iterator(self.filename):
            if line.startswith('#') or line == '':
                continue
            yield Gff3Record.from_string(line)


class _FeatureWriter:
    df: TextIO

    @staticmethod
    def write_iterator(
            iterable: Union[_FeatureIterator, Iterator[Feature]],
            output_filename: str,
            prefix_annotations: Optional[List[str]] = None
    ):
        with _FeatureWriter(output_filename) as writer:
            if prefix_annotations is not None:
                for annotation in prefix_annotations:
                    writer.write_comment(annotation)
            for feature in iterable:
                writer.write_feature(feature)

    def __init__(self, output_filename: str):
        self.output_filename = output_filename
        self.fd = ioctl.get_writer(self.output_filename)

    def write_feature(self, feature: Feature):
        self.fd.write(str(feature) + "\n")

    def write_comment(self, comment: str):
        self.fd.write('#' + comment + "\n")

    def close(self):
        self.fd.close()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return


class GtfWriter(_FeatureWriter):
    pass


class Gff3Writer(_FeatureWriter):
    pass
