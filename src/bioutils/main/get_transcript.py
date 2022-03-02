"""
get_transcript.py -- Filter GTF records by a specific attributes
"""

import argparse
from typing import List

from bioutils.main.sample_transcript import subset_gtf_by_transcript_id
from commonutils.io.tqdm_reader import get_tqdm_line_reader
from commonutils.stdlib_helper.logger_helper import get_logger

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
    for line in get_tqdm_line_reader(args.field_value):
        line = line.strip().strip("\"\'")  # Get rid of quotation marks produced by R
        possible_values.append(line)
    lh.info(f"{len(possible_values)} values loaded")
    subset_gtf_by_transcript_id(
        possible_values=possible_values,
        field_name=args.field_name,
        gtf_filename=args.gtf,
        out_filename=args.out
    )
