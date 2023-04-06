"""
generate_isoform_replicates.py -- Generate Technical Replicates using YASIM V3 API.
"""

__all__ = (
    "create_parser",
    "main"
)

import argparse
from typing import List

import yasim.helper.depth_io
from yasim.helper.frontend import patch_frontend_argument_parser
from yasim.helper import depth


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m yasim generate_isoform_replicates",
        description=__doc__.splitlines()[1]
    )
    parser = patch_frontend_argument_parser(parser, "-d")
    parser.add_argument(
        '-n',
        '--num_replicates',
        required=False,
        help="Number of Replicates to be generated",
        nargs='?',
        type=int,
        action='store',
        default=3
    )
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
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    depth_data = yasim.helper.depth_io.read_depth(args.depth)
    for i in range(args.num_replicates):
        yasim.helper.depth_io.write_depth(depth.generate_depth_replicates_uniform(
            depth_data,
            args.range
        ), f"{args.depth}.{i}", )
