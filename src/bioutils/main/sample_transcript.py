"""
get_transcript.py -- Filter GTF records by a specific attributes
"""

import argparse
import random
from typing import List, Iterable

from bioutils.io.feature import GtfIterator, GtfWriter
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def subset_gtf_by_transcript_id(
        possible_values:Iterable[str],
        field_name:str,
        gtf_filename:str,
        out_filename:str
):
    gi = GtfIterator(gtf_filename)
    final_record_num = 0
    input_record_num = 0
    with GtfWriter(out_filename) as writer:
        for gtf_record in gi:
            input_record_num += 1
            if gtf_record.attribute.get(field_name, None) in possible_values:
                writer.write_feature(gtf_record)
                final_record_num += 1
    lh.info(
        f"{input_record_num} processed with {final_record_num} ({round(final_record_num / input_record_num, 2) * 100}%) records output")


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True, help="Gtf to filter from", nargs='?', type=str, action='store')
    parser.add_argument("--percent", required=False, help="How many percent of transcript to be sampled", nargs='?', type=int,
                        action='store', default=50)
    parser.add_argument("--out", required=True, help="Filtered output", nargs='?', type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gi = GtfIterator(args.gtf)
    transcript_ids = set()
    for gtf_record in gi:
        tmp_tid = gtf_record.attribute.get("transcript_id", None)
        if tmp_tid is not None:
            transcript_ids.add(tmp_tid)
    transcript_ids=random.sample(list(transcript_ids), len(transcript_ids)*args.percent//100)

    subset_gtf_by_transcript_id(
        possible_values=transcript_ids,
        field_name="transcript_id",
        gtf_filename=args.gtf,
        out_filename=args.out
    )

