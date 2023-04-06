"""
generate_barcoded_isoform_replicates.py -- Generate Technical Replicates using YASIM V3 API with barcodes.
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
import os
from typing import List

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from yasim.helper.frontend import patch_frontend_argument_parser
from yasim.helper import depth, depth_io


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m yasim_sc generate_barcoded_isoform_replicates",
        description=__doc__.splitlines()[1]
    )
    parser = patch_frontend_argument_parser(parser, "-d")
    parser = patch_frontend_argument_parser(parser, "-b")
    parser.add_argument(
        '-r',
        '--range',
        required=False,
        help="Range of Generated Data",
        nargs='?',
        type=float,
        action='store',
        default=0.1
    )
    parser.add_argument(
        '-o',
        '--dest_dir_path',
        required=True,
        help="Path of output directory",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    depth_data = depth_io.read_depth(args.depth)
    for barcode in get_tqdm_line_reader(args.barcodes):
        depth_io.write_depth(
            depth.generate_depth_replicates_uniform(
                depth_data,
                args.range
            ),
            os.path.join(args.dest_dir_path, f"{barcode}.tsv"),
            "TRANSCRIPT_ID",
            show_tqdm=False
        )
