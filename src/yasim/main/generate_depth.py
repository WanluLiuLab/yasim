import argparse
from typing import List

from bioutils.datastructure.gene_view import GeneView
from yasim.helper import depth


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output TSV", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--mu', required=False, help="Average depth", nargs='?',
                        type=int, action='store', default=100)
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneView.from_file(args.gtf)
    dge_data = dge_helper.simulate_dge_uniform(
        gv=gv,
        mu=args.mu
    )
    depth.write_dge(dge_data, args.out)
