import argparse
import random
from collections import defaultdict
from typing import Dict, List

from bioutils.datastructure.gene import GeneView
from commonutils import ioctl
from commonutils.tqdm_importer import tqdm

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output TSV", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--max_depth', required=False, help="Max depth", nargs='?',
                        type=int, action='store', default=100)
    parser.add_argument('-l', '--levels', required=False, help="Max number of levels", nargs='?',
                        type=int, action='store', default=100)
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)


def simulate_dge_uniform(
        gv: GeneView,
        output_tsv: str,
        max_depth: int,
        levels: int = 100
) -> Dict[str, int]:
    transcript_names = gv.transcripts.keys()
    depth = {}
    with ioctl.get_writer(output_tsv) as writer:
        writer.write(f"gene_name\tdepth\n")
        for transcript_name in tqdm(iterable=transcript_names, desc="Simulating..."):
            d = int(random.uniform(1, max_depth)) * levels // max_depth * int(max_depth / levels) + 1
            writer.write(f"{transcript_name}\t{d}\n")
            depth[transcript_name] = d
    return depth


def read_depth(input_tsv: str) -> Dict[str, int]:
    retd = {}
    total = ioctl.wc_l(input_tsv)
    with ioctl.get_reader(input_tsv) as reader:
        reader.readline()  # Skip line 1
        for l in tqdm(iterable=reader.readlines(), desc="Reading depth file...", total=total - 1):
            l = l.strip()
            lkv = l.split("\t")
            retd[lkv[0]] = int(lkv[1])
    return retd


def cluster_depth(depth: Dict[str, int]) -> Dict[int, List[str]]:
    retd = defaultdict(lambda: [])
    for k, v in depth.items():
        retd[v].append(k)
    return dict(retd)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneView(args.gtf)
    simulate_dge_uniform(
        gv=gv,
        output_tsv=args.out,
        max_depth=args.max_depth,
        levels=args.levels
    )
