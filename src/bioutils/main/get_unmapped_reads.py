import argparse
from collections import defaultdict
from typing import List

import pysam

from bioutils.datastructure.fastq_view import FastqView

__version__ = 0.1

from bioutils.io.fastq import FastqWriter

from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam', required=True, help="Input SAM/BAM alignment file", nargs='?',
                        type=str, action='store')
    parser.add_argument('-r', '--reads_fq', required=True, help="Input read FASTQ", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output alignment file", nargs='?',
                        type=str, action='store')
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)

    sam = pysam.AlignmentFile(args.sam)
    fqv = FastqView(args.reads_fq)
    unmapped_transcripts_dict = defaultdict(lambda: 0)
    mapped_transcripts_dict = defaultdict(lambda: 0)
    curr_pos = sam.tell()
    full_length = 0
    for _ in sam.fetch(until_eof=True):
        full_length += 1
    sam.seek(curr_pos)
    with FastqWriter(args.out + ".fastq") as unmapped_writer:
        for alignment in tqdm(desc="Aligning unmapped reads...", iterable=sam.fetch(until_eof=True), total=full_length):
            seq_id = alignment.query_name
            if seq_id is None:
                continue
            transcript_name = seq_id.split(":")[0]
            if not alignment.is_unmapped:
                mapped_transcripts_dict[transcript_name] += 1
            else:
                unmapped_transcripts_dict[transcript_name] += 1
                unmapped_writer.write(fqv.get(seq_id))

    with get_writer(args.out + ".stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "SUCCESS_ALN",
            "FAILED_ALN"
        )) + "\n")
        for transcript_id, count in unmapped_transcripts_dict.items():
            stats_writer.write("\t".join((
                transcript_id,
                str(mapped_transcripts_dict[transcript_id]),
                str(count)
            )) + "\n")
