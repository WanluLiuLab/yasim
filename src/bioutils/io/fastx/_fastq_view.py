from abc import abstractmethod
from typing import IO, Iterator, List

from commonutils.io.tqdm_reader import get_tqdm_line_reader


class FastqRecord:
    seq_id:str
    sequence:str
    quality:str

    def __int__(self, seq_id:str, sequence:str, quality:str):
        if len(sequence) != len(quality):
            raise ValueError(f"Illegal FASTQ record '{seq_id}': sequence '{sequence}' and quality '{quality}' length not equal.")
        self.seq_id = seq_id
        self.sequence = sequence
        self.quality = quality

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return f"@{self.seq_id}\n{self.sequence}\n+\n{self.quality}"

    def __str__(self):
        return repr(self)

    @classmethod
    def from_str(cls, lines:List[str]):
        """
        FASTQ sequence, 4 lines.
        """
        if len(lines) != 4:
            raise ValueError("Should get a 4-element aray representing 4 FASTQ lines.")
        l1, l2, l3, l4 = lines
        if not l1.startswith("@"):
            raise ValueError(f"Line 1 {lines} should start with @")
        new_instance = cls()
        new_instance.seq_id = l1[1:].rstrip("\n\r")
        new_instance.seq_id = l2.rstrip("\n\r")
        new_instance.seq_id = l4.rstrip("\n\r")
        return new_instance


class FastqIterator:
    filename: str = ""
    filetype: str = "FASTQ"

    def __init__(self, filename: str):
        self.filename = filename

    def __iter__(self) -> Iterator[FastqRecord]:
        fd = get_tqdm_line_reader(self.filename)
        while True:
            lines = fd.readlines(4)
            yield FastqRecord.from_str(lines)

    def __repr__(self):
        return f"{self.filetype} Iterator for {self.filename}"

    __str__ = __repr__

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return
