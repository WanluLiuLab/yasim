import argparse
from typing import List

from yasim.helper.tcr import create_tcr_cache


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--tcr_genelist_path',
        required=True,
        help="TCR Gene List JSON",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-f',
        '--fasta',
        required=True,
        help="Path to reference genome FASTA",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-g',
        '--gtf',
        required=True,
        help="Path to reference genome GTF",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--tcr_aa_table_path',
        required=True,
        help="TCR AA Table Path",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output TCR Cache",
        nargs='?',
        type=str,
        action='store'
    )
    return parser.parse_args(args)


def main(args: List[str]) -> int:
    args = _parse_args(args)
    basename = "."
    create_tcr_cache(
        ref_fa_path=args.fasta,
        ref_gtf_path=args.gtf,
        tcr_genelist_path=args.tcr_genelist_path,
        tcr_aa_table_path=args.tcr_aa_table_path,
        tcr_cache_path=args.out
    )
    return 0
