import argparse
from typing import List

from yasim.helper import depth


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depth', required=True, help="Depth TSV", nargs='?',
                        type=str, action='store')
    parser.add_argument('-n', '--num_replicates', required=False, help="Number of Replicates", nargs='?',
                        type=int, action='store', default=3)
    parser.add_argument('-r', '--range', required=False, help="Range of Generated Data", nargs='?',
                        type=float, action='store', default=0.001)
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    depth_data = depth.read_depth(args.depth)
    for i in range(args.num_replicates):
        depth.write_depth(
            depth.generate_depth_replicates_uniform(
                depth_data,
                args.range
            ),
            f"{depth_data}.{i}"
        )
