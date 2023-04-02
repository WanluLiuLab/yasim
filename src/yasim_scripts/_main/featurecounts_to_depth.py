import argparse
from typing import List

import pandas as pd

from yasim.helper.depth_io import write_depth


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--fc_tsv',
        required=True,
        help="featureCounts output TSV",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output TSV",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-f',
        '--feature_name',
        required=False,
        help="Name of output feature",
        choices=("GENE_ID", "TRANSCRIPT_ID"),
        default="TRANSCRIPT_ID",
        nargs='?',
        type=str,
        action='store'
    )
    return parser.parse_args(args)


def main(args: List[str]) -> int:
    depth = {}
    args = _parse_args(args)
    fc = pd.read_table(
        args.fc_tsv,
        sep="\t",
        comment="#"
    )
    for datum in fc.itertuples(index=False):
        real_depth = datum[-1] # round(datum[-1] / datum[5], 2)
        if real_depth < 0.01:
            continue
        depth[datum[0]] = real_depth
    write_depth(depth, args.out, args.feature_name)
    return 0
