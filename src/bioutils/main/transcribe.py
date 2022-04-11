import argparse
from typing import List

from bioutils.datastructure.fasta_view import FastaViewFactory
from bioutils.datastructure.gene_view import GeneViewFactory
from bioutils.datastructure.gv_helper import transcribe
from commonutils.stdlib_helper.logger_helper import get_logger

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Name of Output FASTA", nargs='?',
                        type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    fv = FastaViewFactory(args.fasta)
    transcribe(gv, args.out, fv)
