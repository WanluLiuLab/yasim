from typing import Iterator, TextIO, Union

from bioutils.typing.fastq import FastqRecord
from bioutils.typing.misc import BaseIterator
from commonutils.io.safe_io import get_writer, get_reader
from commonutils.io.tqdm_reader import get_tqdm_line_reader


class FastqIterator(BaseIterator):
    filetype: str = "FASTQ"
    record_type = FastqRecord

    def __init__(self, filename: str, show_tqdm: bool = True):
        super().__init__(filename)
        if show_tqdm:
            self.fd = get_tqdm_line_reader(self.filename)
        else:
            self.fd = get_reader(self.filename)

    def __iter__(self) -> Iterator[FastqRecord]:
        while True:
            lines = [self.fd.readline(-1) for _ in range(4)]
            if '' in lines:
                break
            yield FastqRecord.from_str(lines)


class FastqWriter:
    df: TextIO

    @staticmethod
    def write_iterator(
            iterable: Union[FastqIterator, Iterator[FastqRecord]],
            output_filename: str
    ):
        with FastqWriter(output_filename) as writer:
            for fastq_record in iterable:
                writer.write(fastq_record)

    def __init__(self, output_filename: str):
        self.output_filename = output_filename
        self.fd = get_writer(self.output_filename)

    def write(self, fastq_record: FastqRecord):
        self.fd.write(str(fastq_record) + "\n")

    def close(self):
        self.fd.close()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return

    def tell(self) -> int:
        try:
            return self.fd.tell()
        except OSError:
            return -1

    def __repr__(self) -> str:
        return f"FASTQ writer of {self.output_filename} at {self.tell()}"

    __str__ = __repr__
