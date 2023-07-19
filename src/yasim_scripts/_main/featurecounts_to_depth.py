"""
featurecounts_to_depth.py -- Convert FeatureCounts Output for NGS and TGS to YASIM input depth.

.. versionadded:: 3.1.5
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse

import pandas as pd

from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List
from yasim.helper.depth_io import write_depth


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_scripts featurecounts_to_depth",
        description=__doc__.splitlines()[1]
    )
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help="Path to featureCounts/Salmon output TSV",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--software',
        required=False,
        help="name of quantification software",
        nargs='?',
        choices=("featureCounts", "Salmon"),
        default="featureCounts",
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Path to output TSV",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-f',
        '--feature_name',
        required=False,
        help="Name of output feature",
        choices=("GENE_ID", "TRANSCRIPT_ID"),
        default="TRANSCRIPT_ID",
        nargs='?',
        type=str,
        action='store'
    )
    ngs_tgs_meg = parser.add_mutually_exclusive_group()
    ngs_tgs_meg.add_argument(
        '--read_length',
        required=False,
        help="[For NGS Only] Read length",
        nargs='?',
        type=int,
        default=None,
        action='store'
    )
    ngs_tgs_meg.add_argument(
        '--read_completeness',
        required=False,
        help="[For TGS Only] Mean read completeness",
        nargs='?',
        type=float,
        default=None,
        action='store'
    )
    return parser


def main(args: List[str]) -> int:
    depth = {}
    args = create_parser().parse_args(args)
    fc = pd.read_table(
        args.input,
        sep="\t",
        comment="#"
    )
    if args.read_length is not None:
        if args.software == "featureCounts":
            def get_read_depth(_datum):
                return _datum[-1] * args.read_length / _datum[5]
        else:
            def get_read_depth(_datum):
                return _datum.NumReads * args.read_length / _datum.Length
    elif args.read_completeness is not None:
        def get_read_depth(_datum):
            return _datum[-1] * args.read_completeness
    else:
        def get_read_depth(_datum):
            return _datum[-1]

    for datum in fc.itertuples(index=False):
        real_depth = get_read_depth(datum)
        print(real_depth, datum)
        if real_depth < 0.01:
            continue
        depth[datum[0]] = real_depth
    write_depth(depth, args.out, args.feature_name)
    return 0
