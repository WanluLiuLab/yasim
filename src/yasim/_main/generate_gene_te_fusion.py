"""
insert_transposon.py -- Generate 

.. versionadded:: 3.2.0
"""

__all__ = ("main", "create_parser")

import argparse
import random

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_tree import GeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List
from yasim.helper.frontend import patch_frontend_argument_parser
from yasim.helper.translation_instruction import TranslationInstruction
from yasim.helper.transposon import TransposonDatabase

rdg = random.SystemRandom()


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim insert_transposon",
        description=__doc__.splitlines()[1],
    )
    parser = patch_frontend_argument_parser(parser, "-g")
    parser = patch_frontend_argument_parser(parser, "-f")
    parser.add_argument(
        "--tedb",
        type=str,
        help="Path to TE database index.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        help="Path to base path where result file should be written to.",
        required=True,
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
    argv = create_parser().parse_args(args)
    ti = TranslationInstruction.generate(
        n=1000,
        mu=argv.mu,
        fav=FastaViewFactory(argv.fasta),
        gt=GeneTree.from_gtf_file(argv.gtf, gene_implementation=DumbGene),
        tedb=TransposonDatabase.load(argv.tedb, with_tqdm=True),
    )
    ti.to_fasta(argv.out + ".fa")
    ti.to_depth(argv.out + ".depth.tsv")
    ti.to_json(argv.out + ".json")
