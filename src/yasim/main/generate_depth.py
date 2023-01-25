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
                        type=float, action='store', default=100)
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    dge_data = depth.simulate_dge_gmm(
        gv=gv,
        mu=args.mu
    )
    depth.write_dge(dge_data, args.out)
