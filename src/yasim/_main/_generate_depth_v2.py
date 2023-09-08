"""
generate_depth_v2.py -- Generate Isoform-Level Depth using YASIM V2 API.

.. versionadded:: 0.2.7
    Named as generate_depth

.. versionchanged:: 3.0.0
    Renamed to generate_depth_v2

.. versiondeprecated:: 3.1.5
    Use V3API instead.
"""

__all__ = ("main", "create_parser")

import argparse

from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from yasim.helper import depth, depth_io
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim generate_depth_v2",
        description=__doc__.splitlines()[1],
    )
    parser = patch_frontend_argument_parser(parser, "-g")
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to output Isoform-Level depth TSV. Can be compressed.",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-d",
        "--mu",
        required=False,
        help="Average depth",
        nargs="?",
        type=float,
        action="store",
        default=100,
    )
    return parser


def main(args: List[str]):
    _lh.warning(
        "DEPRECATION WARNING: "
        "The V2 API had been deprecated. Use `generate_gene_depth` and `generate_isoform_depth` instead"
    )
    args = create_parser().parse_args(args)
    gv = DiploidGeneTree.from_gtf_file(args.gtf, gene_implementation=DumbGene)
    dge_data = depth.simulate_depth_gmm_v2(isoform_names=gv.transcript_ids, mu=args.mu)
    depth_io.write_depth(dge_data, args.out, feature_name="TRANSCRIPT_ID")
