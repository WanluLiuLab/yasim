from collections import defaultdict

from bioutils.algorithm.alignment import smith_waterman_backtrack
from bioutils.datastructure.fasta_view import FastaView
from bioutils.datastructure.fastq_view import FastqView

import argparse
from typing import List

__version__ = 0.1

from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam', required=True, help="Input SAM/BAM alignment file", nargs='?',
                        type=str, action='store')
    parser.add_argument('-f', '--fasta', required=True, help="Input transcriptom FASTA", nargs='?',
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
    fav = FastaView(args.fasta)
    unmapped_transcripts_dict = defaultdict(lambda: 0)
    curr_pos = sam.tell()
    full_length = 0
    for _ in sam.fetch(until_eof=True):
        full_length += 1
    sam.seek(curr_pos)
    with get_writer(args.out) as aln_writer:
        for alignment in tqdm(iterable=sam.fetch(until_eof=True), total=full_length):
            if not alignment.is_unmapped:
                pass
            seq_id = alignment.query_name
            if seq_id is None:
                continue
            fastq_record = fqv.get(seq_id)
            transcript_name = seq_id.split(":")[0]
            sw_backtrack = smith_waterman_backtrack(
                seq1=fastq_record.sequence,
                seq2=fav.sequence(seq_id),
                alignment_title=seq_id
            )
            unmapped_transcripts_dict[transcript_name] += 1
        for backtrack in sw_backtrack:
            aln_writer.write(backtrack + "/n")
    with get_writer(args.out + ".stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "ALN",
            "FAILED_ALN"
        )))
        for transcript_id, count in unmapped_transcripts_dict:
            pass
