




from typing import Iterator

from bioutils.typing.fastq import FastqRecord
from bioutils.typing.misc import BaseIterator
from commonutils.io.tqdm_reader import get_tqdm_line_reader


class FastqIterator(BaseIterator):
    filetype: str = "FASTQ"
    record_type = FastqRecord

    def __iter__(self) -> Iterator[FastqRecord]:
        fd = get_tqdm_line_reader(self.filename)
        while True:
            lines = fd.readlines(4)
            yield FastqRecord.from_str(lines)
