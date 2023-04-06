"""
generate_as_events.py -- Generate Alternative Splicing Events from Reference genome using YASIM V3 API.
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
from typing import List

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.as_events import ASManipulator
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m yasim generate_as_events", description=__doc__.splitlines()[1])
    parser = patch_frontend_argument_parser(parser, "-g")
    parser = patch_frontend_argument_parser(parser, "-f")
    parser.add_argument(
        '-c',
        '--complexity',
        required=True,
        help="Genome Complexity Index, should be an integer between 1 and 9.",
        nargs='?',
        type=int,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Path to output genome annotation with Ground-Truth AS Events in GTF.",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    _lh.info(f"Loaded {gv.number_of_genes} genes with {gv.number_of_transcripts} transcript")
    asm = ASManipulator(gv=gv)
    asm.run("ce", args.complexity)
    asm.to_file(args.out)
