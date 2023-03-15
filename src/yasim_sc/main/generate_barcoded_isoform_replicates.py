import argparse
import os
from typing import List

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from yasim.helper import depth


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depth', required=True, help="Depth TSV", nargs='?',
                        type=str, action='store')
    parser.add_argument('-b', '--barcodes', required=True, help="Barcode TXT", nargs='?',
                        type=str, action='store')
    parser.add_argument('-r', '--range', required=False, help="Range of Generated Data", nargs='?',
                        type=float, action='store', default=0.1)
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    depth_data = depth.read_depth(args.depth)
    for barcode in get_tqdm_line_reader(args.barcodes):
        depth.write_depth(
            depth.generate_depth_replicates_uniform(
                depth_data,
                args.range
            ),
            os.path.join(f"{args.depth}.d", f"{barcode}.tsv")
        )
