import argparse
from typing import List

from labw_utils.bioutils.datastructure.gene_view import GeneViewFactory
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

from yasim.helper.as_events import ASManipulator

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Reference genome, in GTF format", nargs='?',
                        type=str, action='store')
    parser.add_argument(
        '-c', '--complexity',
        required=True,
        help="Genome Complexity, should be an integer between 1 and 9",
        nargs='?',
        type=int,
        action='store'
    )
    parser.add_argument('-o', '--out', required=True, help="Output GTF", nargs='?',
                        type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    logger.info(f"Loaded {gv.number_of_genes} genes with {gv.number_of_transcripts} transcript")
    asm = ASManipulator(gv=gv)
    asm.run("ce", args.complexity)  # TODO: add more organisms
    asm.to_file(args.out)
