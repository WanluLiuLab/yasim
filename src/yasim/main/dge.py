import argparse
import random
from typing import Dict, List

from bioutils.datastructure.gene_view import GeneView
from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer, get_reader

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
    transcript_ids = gv.transcripts.keys()
    depth = {}
    with get_writer(output_tsv) as writer:
        writer.write(f"TRANSCRIPT_ID\tDEPTH\n")
        for transcript_id in tqdm(iterable=transcript_ids, desc="Simulating..."):
            d = int(random.uniform(1, max_depth)) * levels // max_depth * int(max_depth / levels) + 1
            writer.write(f"{transcript_id}\t{d}\n")
            depth[transcript_id] = d
    return depth


def read_depth(input_tsv: str) -> Dict[str, int]:
    retd = {}
    total = shell_utils.wc_l(input_tsv)
    with get_reader(input_tsv) as reader:
        reader.readline()  # Skip line 1
        for line in tqdm(iterable=reader.readlines(), desc="Reading depth file...", total=total - 1):
            line = line.strip()
            lkv = line.split("\t")
            retd[lkv[0]] = int(lkv[1])
    return retd


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneView.from_file(args.gtf)
    simulate_dge_uniform(
        gv=gv,
        output_tsv=args.out,
        max_depth=args.max_depth,
        levels=args.levels
    )
