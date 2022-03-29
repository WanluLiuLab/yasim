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