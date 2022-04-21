"""
gtf_head -- Get first N records in GTF
"""

import argparse
import math
from typing import List

from bioutils.datastructure.gene_view import GeneViewFactory
from bioutils.io.feature import GtfWriter
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF", nargs='?', type=str, action='store')
    parser.add_argument("--gene", required=False, help="Get first N genes", nargs='?', type=int, action='store',
                        default=math.inf)
    parser.add_argument("--transcript", required=False, help="Get first N transcripts", nargs='?', type=int,
                        action='store', default=math.inf)
    parser.add_argument("--record", required=False, help="Get first N records", nargs='?', type=int, action='store',
                        default=math.inf)

    parser.add_argument("-o", "--out", required=True, help="Output GTF", nargs='?', type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf, not_save_index=True)
    gi = gv.get_iterator()
    n_of_record = 0
    n_of_gene = 0
    n_of_transcript = 0
    with GtfWriter(args.out) as gw:
        for feature in gi:
            n_of_record += 1
            if feature.feature == "gene":
                n_of_gene += 1
            elif feature.feature in ("transcript", "mRNA"):
                n_of_transcript += 1
            if n_of_gene > args.gene or n_of_transcript > args.transcript or n_of_record > args.record:
                break
            gw.write_feature(feature)
