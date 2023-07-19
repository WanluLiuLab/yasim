"""
generate_as_events.py -- Generate Alternative Splicing Events from Reference genome using YASIM V3 API.

.. versionadded:: 3.1.5
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from yasim.helper.as_events import ASManipulator
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(prog="python -m yasim generate_as_events",
                                                  description=__doc__.splitlines()[1])
    parser = patch_frontend_argument_parser(parser, "-g")
    parser = patch_frontend_argument_parser(parser, "-f")
    parser.add_argument(
        '-c',
        '--complexity',
        required=True,
        help="Transcriptome Complexity Index, should be an integer between 1 and 9.",
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
    _lh.info(
        "GEN AS EVENTS: Loaded %d genes with %d transcript",
        gv.number_of_genes,
        gv.number_of_transcripts
    )
    asm = ASManipulator(gv=gv)
    asm.run("ce", args.complexity)
    asm.to_file(args.out)
