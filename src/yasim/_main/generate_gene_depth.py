"""
generate_gene_depth.py -- Generate Gene-Level Depth using YASIM V3 API.

.. versionadded:: 3.1.5
"""

__all__ = ("main", "create_parser")

import argparse
from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene

from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from yasim.helper import depth, depth_io
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim generate_gene_depth",
        description=__doc__.splitlines()[1],
    )
    parser = patch_frontend_argument_parser(parser, "-g")
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to output Gene-Level depth TSV. Can be compressed.",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-d",
        "--mu",
        required=False,
        help="Average depth.",
        nargs="?",
        type=float,
        action="store",
        default=depth.DEFAULT_MU,
    )
    parser = patch_frontend_argument_parser(parser, "--low_cutoff")
    parser = patch_frontend_argument_parser(parser, "--high_cutoff_ratio")
    return parser


def main(args: List[str]) -> int:
    args = create_parser().parse_args(args)
    gv = DiploidGeneTree.from_gtf_file(args.gtf, gene_implementation=DumbGene)
    try:
        dge_data = depth.simulate_gene_level_depth_gmm(
            gene_names=gv.gene_ids,
            mu=args.mu,
            low_cutoff=args.low_cutoff,
            high_cutoff_ratio=args.high_cutoff_ratio,
        )
    except depth.GenerationFailureException:
        _lh.error("GEN GENE DEPTH: Generation failed! You may try again or use larger `mu`.")
        return 1
    depth_io.write_depth(dge_data, args.out, feature_name="GENE_ID")
