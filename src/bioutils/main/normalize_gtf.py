"""
normalize_gtf -- Performs GTF normalization, etc.
"""

import argparse
from typing import List

from bioutils.datastructure.gene_view import GeneViewFactory
from bioutils.datastructure.gene_view_proxy import DEFAULT_SORT_EXON_EXON_STRAND_POLICY, \
    VALID_SORT_EXON_EXON_STRAND_POLICY
from bioutils.io.feature import GtfIterator, GtfWriter
from bioutils.typing.feature import VAILD_GTF_QUOTE_OPTONS
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF", nargs='?', type=str, action='store')
    parser.add_argument("--three_tier",
                        help="Whether to parse the GTF into Gene-Transcript-Exon Three-Tier Structure. Other features will be discarded, and missing genes & transcripts will be added by maximun span length",
                        action='store_true')
    parser.add_argument("--quote", required=False, help="Whether to add quotes in alternative field of output GTF",
                        nargs='?', type=str, action='store', choices=VAILD_GTF_QUOTE_OPTONS, default="all")
    parser.add_argument(
        "--sort_exon_exon_number_policy", required=False, help="How you would rank exon_number? On Reference: unstranded. Sorted by bedtools: stranded",
        nargs='?', type=str, action='store', choices=VALID_SORT_EXON_EXON_STRAND_POLICY, default=DEFAULT_SORT_EXON_EXON_STRAND_POLICY
    )
    parser.add_argument("--out", required=True, help="Output GTF", nargs='?', type=str, action='store')
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    if args.three_tier:
        gv = GeneViewFactory.from_file(args.gtf, not_save_index=True)
        gv.standardize(
            sort_exon_exon_number_policy=args.sort_exon_exon_number_policy
        )
        gi = gv.get_iterator()
    else:
        gi = GtfIterator(args.gtf)
    with GtfWriter(args.out) as gw:
        for feature in gi:
            gw.write_feature(feature, quote=args.quote)
