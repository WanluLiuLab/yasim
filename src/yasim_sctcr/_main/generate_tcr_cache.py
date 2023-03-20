import argparse
from typing import List

from yasim.helper.tcr import create_tcr_cache
from yasim_sctcr._main import get_sample_data_path


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--tcr_genelist_path',
        required=True,
        help=f"TCR Gene List JSON. The IMGT version for human is {get_sample_data_path('tcr_gene_list.min.json')}",
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
        help=f"TCR AA Table Path. The IMGT version for human is {get_sample_data_path('tcr_aa_table.min.json')}",
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
    create_tcr_cache(
        ref_fa_path=args.fasta,
        ref_gtf_path=args.gtf,
        tcr_genelist_path=args.tcr_genelist_path,
        tcr_aa_table_path=args.tcr_aa_table_path,
        tcr_cache_path=args.out
    )
    return 0
