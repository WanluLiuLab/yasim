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
            gene.number_of_transcripts,
            gene_level_depth[gene.gene_id]
        )
        for i, transcript in enumerate(gene.iter_transcripts()):
            transcript_level_depth[transcript.transcript_id] = this_transcript_level_depth[i]
    depth.write_depth(transcript_level_depth, args.out)
