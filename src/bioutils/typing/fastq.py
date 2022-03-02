"""
fastq.py -- An Im-Memory FastQ Record.

>>> fastq_seq = ["@A\\n", "AGCT\\n", "+\\n", "0000\\n"]
>>> fastq_record = FastqRecord.from_str(fastq_seq)
>>> fastq_record.seq_id
'A'
>>> fastq_record.sequence
'AGCT'
>>> fastq_record.quality
'0000'
>>> str(fastq_record) == '@A\\nAGCT\\n+\\n0000'
True
>>> fastq_record = FastqRecord.from_single_str('@A\\nAGCT\\n+\\n0000')
>>> str(fastq_record) == '@A\\nAGCT\\n+\\n0000'
True
"""

from typing import List


class FastqRecord:
    """
    A naive in-memory FASTQ record.
    """

    seq_id: str
    """
    Sequence ID.
    """

    sequence: str
    """
    The sequence.
    """

    quality: str
    """
    The corresponding quality, whose length should be equal to ``sequence``.
    """

    def __int__(self, seq_id: str, sequence: str, quality: str):
        if len(sequence) != len(quality):
            raise ValueError(
                f"Illegal FASTQ record '{seq_id}': sequence '{sequence}' and quality '{quality}' length not equal.")
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
    def from_str(cls, lines: List[str]):
        """
        Generate from FASTQ sequence, 4 lines.

        This method is set to generate record from arbitrary :py:mod:`typing.TextIO` readers.
        """
        if len(lines) != 4:
            raise ValueError("Should get a 4-element aray representing 4 FASTQ lines.")
        l1, l2, l3, l4 = lines
        if not l1.startswith("@"):
            raise ValueError(f"Line 1 {lines} should start with @")
        new_instance = cls()
        new_instance.seq_id = l1[1:].rstrip("\n\r")
        new_instance.sequence = l2.rstrip("\n\r")
        new_instance.quality = l4.rstrip("\n\r")
        return new_instance

    @classmethod
    def from_single_str(cls, input_str: str):
        """
        Generate from FASTQ sequence, 1 line.
        """
        return cls.from_str(input_str.splitlines())
