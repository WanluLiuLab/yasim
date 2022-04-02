import argparse
import random
from typing import List

from bioutils.datastructure.fasta_view import FastaViewFactory, FastaViewType
from bioutils.datastructure.gene_view import GeneViewFactory, GeneViewType
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Reference genome, in GTF format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output GTF", nargs='?',
                        type=str, action='store')
    return parser.parse_args(args)


def sample_exon(
        gv: GeneViewType,
        output_gtf_filename: str,
        fasta_handler: FastaViewType
) -> GeneViewType:
    transcript_name_to_del = []
    for k, v in tqdm(iterable=gv.transcripts.items(), desc="Sampling Exons..."):
        indices = random.sample(range(len(v.exons)), int(len(v.exons) * 0.75))
        v.exons = [v.exons[i] for i in sorted(indices)]
        if len(v.cdna_sequence(sequence_func=fasta_handler.sequence)) >= 250:
            pass
        else:
            transcript_name_to_del.append(k)
    for transcript_name in tqdm(iterable=transcript_name_to_del, desc="Deleting unwanted transcripts"):
        gv.del_transcript(transcript_name)
    logger.info(f"Remaining {len(gv.genes)} genes with {len(gv.transcripts)} transcript")
    gv.standardize()
    gv.to_file(output_gtf_filename)
    return gv


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    logger.info(f"Loaded {len(gv.genes)} genes with {len(gv.transcripts)} transcript")
    fv = FastaViewFactory(args.fasta)
    sample_exon(gv=gv, output_gtf_filename=args.out, fasta_handler=fv)
