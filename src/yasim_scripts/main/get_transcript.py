"""
get_transcript.py -- Filter GTF records by a specific attributes
"""

import argparse
from typing import List

from bioutils.io.feature import GtfIterator, GtfWriter
from commonutils.io.safe_io import get_reader
from commonutils.logger import get_logger

lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True, help="Gtf to filter from", nargs='?', type=str, action='store')
    parser.add_argument("--field_name", required=False, help="Attribute to be filtered", nargs='?', type=str,
                        action='store', default="transcript_id")
    parser.add_argument("--field_value", required=True,
                        help="Filename that contains legal values, once per line, can be quoted", nargs='?', type=str,
                        action='store')
    parser.add_argument("--out", required=True, help="Filtered output", nargs='?', type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    possible_values = []
    with get_reader(args.field_value) as reader:
        while True:
            line = reader.readline()
            if not line:
                break
            line = line.strip().strip("\"\'")  # Get rid of quotation marks produced by R
            possible_values.append(line)
    lh.info(f"{len(possible_values)} values loaded")
    gi = GtfIterator(args.gtf)
    final_record_num = 0
    input_record_num = 0
    with GtfWriter(args.out) as writer:
        for gtf_record in gi:
            input_record_num += 1
            if gtf_record.attribute.get(args.field_name, None) in possible_values:
                writer.write_feature(gtf_record)
                final_record_num += 1
    lh.info(
        f"{input_record_num} processed with {final_record_num} ({round(final_record_num / input_record_num, 2) * 100}%) records output")
