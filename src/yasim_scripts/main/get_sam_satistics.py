"""A General-Purposed SAM quality control tool"""

import argparse
from statistics import mean
from typing import List, Any

import pysam

from bioutils.algorithm.sequence import get_gc_percent
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    sam = pysam.AlignmentFile(args.sam)
    num_mapped_reads = 0

    with get_writer(args.out) as writer:
        writer.write(
            "\t".join((
                "LEN",
                "GC",
                "IS_MAPPED",
                "MAPQ",
                "ALNQ"
            )) + "\n"
        )
        curr_pos = sam.tell()
        full_length = 0
        for _ in sam.fetch(until_eof=True):
            full_length += 1
        sam.seek(curr_pos)

        for alignment in tqdm(iterable=sam.fetch(until_eof=True), total=full_length):
            if alignment.is_unmapped:
                ismapped = 0
            else:
                ismapped = 1
            num_mapped_reads += ismapped
            if alignment.query_length is None or alignment.query_qualities is None or alignment.query_sequence is None:
                continue
            read_lenth = alignment.query_length
            alnq = mean(alignment.query_qualities)
            writer.write(
                "\t".join((
                    str(read_lenth),
                    str(round(get_gc_percent(alignment.query_sequence) * 100, 2)),
                    str(ismapped),
                    str(round(alignment.mapping_quality, 2)),
                    str(round(alnq, 2))
                )) + "\n"
            )
    with get_writer(args.out + ".summary") as writer:
        def k_v_write(k: str, v: Any):
            rv = repr(v)
            if len(rv) >= 2 and rv[0] == '\'' and rv[-1] == '\'':
                rv = rv[1:-1]
            writer.write("\t".join((k, rv)) + "\n")

        k_v_write("KEY", "VALUE")
        k_v_write("N_READS", full_length)
        k_v_write("MAPPING_RATE", round(num_mapped_reads / full_length * 100, 3))
