"""
generate_isoform_depth.py -- Generate Isoform-Level Depth using YASIM V3 API.
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
from labw_utils.typing_importer import List

import yasim.helper.depth_io
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper import depth
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(prog="python -m yasim generate_isoform_depth", description=__doc__.splitlines()[1])
    parser = patch_frontend_argument_parser(parser, "-g")
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Path to output Isoform-Level Depth TSV. Can be compressed.",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        help="Path to input Gene-Level Depth TSV. Can be compressed.",
        nargs='?',
        type=str,
        action='store'
    )
    parser = patch_frontend_argument_parser(parser, "--low_cutoff")
    parser.add_argument(
        '--alpha',
        required=False,
        help="Zipf's Coefficient, larger for larger differences",
        nargs='?',
        type=int,
        action='store',
        default=4
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    gene_level_depth = yasim.helper.depth_io.read_depth(args.depth)
    transcript_level_depth = {}
    for gene in gv.iter_genes():
        if gene.gene_id not in gene_level_depth:
            _lh.warning("Gene %s defined in GTF but not gene-level depth", gene.gene_id)
        if gene_level_depth[gene.gene_id] == 0:
            for transcript in gene.iter_transcripts():
                transcript_level_depth[transcript.transcript_id] = 0
            continue
        try:
            this_transcript_level_depth = depth.simulate_isoform_variance_inside_a_gene(
                n=gene.number_of_transcripts,
                mu=gene_level_depth[gene.gene_id],
                low_cutoff=args.low_cutoff,
                alpha=args.alpha
            )
        except depth.GenerationFailureException:
            _lh.error("Generation failed for gene %s -- SKIPPED", gene.gene_id)
            continue
        for i, transcript in enumerate(gene.iter_transcripts()):
            transcript_level_depth[transcript.transcript_id] = this_transcript_level_depth[i]
    yasim.helper.depth_io.write_depth(transcript_level_depth, args.out, "TRANSCRIPT_ID")
