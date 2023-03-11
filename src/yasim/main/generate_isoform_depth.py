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
    parser.add_argument('-d', '--depth', required=False, help="Input Gene-Level Expression Abundance", nargs='?',
                        type=str, action='store', default=100)
    parser.add_argument('--low_cutoff', required=False, help="Depth lower than this value would be 0.", nargs='?',
                        type=float, action='store', default=0.01)
    parser.add_argument('--alpha', required=False, help="Zipf's Coefficient, larger for larger differences", nargs='?',
                        type=int, action='store', default=4)
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    gene_level_depth = depth.read_depth(args.depth)
    transcript_level_depth = {}
    for gene in gv.iter_genes():
        if gene_level_depth[gene.gene_id] == 0:
            for transcript in gene.iter_transcripts():
                transcript_level_depth[transcript.transcript_id] = 0
            continue
        this_transcript_level_depth = depth.simulate_isoform_variance_inside_a_gene(
            n=gene.number_of_transcripts,
            mu=gene_level_depth[gene.gene_id],
            low_cutoff=args.low_cutoff,
            alpha=args.alpha
        )
        for i, transcript in enumerate(gene.iter_transcripts()):
            transcript_level_depth[transcript.transcript_id] = this_transcript_level_depth[i]
    depth.write_depth(transcript_level_depth, args.out)
