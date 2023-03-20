__all__ = (
    "main",
)

import argparse
from typing import List

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from yasim.helper.depth_io import write_depth


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--barcode_path',
        required=True,
        help="Barcode TXT path generated by `yasim_sc generate_barcode`",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output basename",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        help="Simulated depth",
        nargs='?',
        type=int,
        action='store'
    )
    return parser.parse_args(args)


def main(args: List[str]) -> int:
    args = _parse_args(args)
    depth_db = {}
    for barcode in get_tqdm_line_reader(args.barcode_path):
        depth_db[f"{barcode}:A"] = args.depth
        depth_db[f"{barcode}:B"] = args.depth

    write_depth(depth_db, args.out, "TRANSCRIPT_ID")
    return 0
