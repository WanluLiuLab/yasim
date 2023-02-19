import argparse
from typing import List

from labw_utils.bioutils.datastructure.gene_view import GeneViewFactory

from yasim.helper import depth


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output TSV", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--mu', required=False, help="Average depth", nargs='?',
                        type=float, action='store', default=20)
    parser.add_argument('--low_cutoff', required=False, help="Depth lower than this value would be 0.", nargs='?',
                        type=float, action='store', default=0.01)
    parser.add_argument('--high_cutoff_ratio', required=False,
                        help="Depth higher than mu * high_cutoff_ratio would be 0.", nargs='?',
                        type=float, action='store', default=10)
    return parser.parse_args(args)

def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    dge_data = depth.simulate_gene_level_depth_gmm(
        gv=gv,
        mu=args.mu,
        low_cutoff=args.low_cutoff,
        high_cutoff_ratio=args.high_cutoff_ratio
    )
    depth.write_depth(dge_data, args.out, feature_name="GENE_ID")
