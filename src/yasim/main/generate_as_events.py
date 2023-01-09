import argparse
from typing import List

from yasim.helper.as_events import ASManipulator

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory, FastaViewType
from labw_utils.bioutils.datastructure.gene_view import GeneViewFactory, GeneViewType
from labw_utils.bioutils.datastructure.gv_feature_proxy import Transcript
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger


logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Reference genome, in GTF format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output GTF", nargs='?',
                        type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    logger.info(f"Loaded {gv.number_of_genes} genes with {gv.number_of_transcripts} transcript")
    asm = ASManipulator(gv=gv)
    asm.run("ce")  # TODO: add more organisms
    asm.to_file(args.out)
