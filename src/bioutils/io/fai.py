from typing import Dict, Iterator

from bioutils.typing.fai import FastaIndexEntry
from commonutils.io import determine_file_line_endings
from commonutils.io.safe_io import get_writer
from commonutils.io.tqdm_reader import get_tqdm_line_reader

FAI_INDEX_TYPE = Dict[str, FastaIndexEntry]


def get_fai_iterator_from_fasta(filename: str) -> Iterator[FastaIndexEntry]:
    """
    Get in-memory FAI index for Fasta

    Do not use this feature on full headers.
    """
    name = ""
    offset = 0
    line_len = 0
    line_blen = 0
    length = 0

    def generate_fai_record():
        new_fai_index = FastaIndexEntry()
        new_fai_index.name = name
        new_fai_index.line_len = line_len
        new_fai_index.line_blen = line_blen
        new_fai_index.offset = offset
        new_fai_index.length = length
        return new_fai_index

    with get_tqdm_line_reader(filename, newline=determine_file_line_endings(filename)) as reader:
        while True:
            line = reader.readline()
            if not line:
                break
            if line == '':
                continue
            if line[0] == '>':  # FASTA header
                if name != '':
                    yield generate_fai_record()
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
            yield generate_fai_record()


def create_fai_from_fasta(
        filename: str,
        index_filename: str
) -> FAI_INDEX_TYPE:
    """
    Create an FAI index for Fasta

    Do not use this feature on full headers.
    """
    return_fai: FAI_INDEX_TYPE = {}
    for fai_record in get_fai_iterator_from_fasta(filename):
        return_fai[fai_record.name] = fai_record

    with get_writer(index_filename) as writer:
        for fai_record in return_fai.values():
            writer.write(str(fai_record) + "\n")
    return return_fai
