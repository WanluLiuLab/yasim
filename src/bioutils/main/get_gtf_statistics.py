"""Get statistics about GTF files that can be parsed into a Gene-Transcript-Exon Three-Tier Structure"""

import argparse
from typing import List

from bioutils.datastructure.gene_view import GeneViewFactory, GeneViewType
from bioutils.io.feature import GtfWriter
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True)
    parser.add_argument("-o", "--out", help="Output prefix", required=False, default=None)
    parser.add_argument(
        "--get_overlapping_transcripts_from_different_gene",
        help="",
        action='store_true'
    )
    return parser.parse_args(args)


def get_overlapping_transcripts_from_different_gene(
        out_basename:str,
        gv:GeneViewType
):
    transcripts = list(gv.iter_transcripts())
    with GtfWriter(f"{out_basename}.overlapping_transcript.gtf") as writer:
        for t_i in tqdm(desc="Iterating over transcripts...", iterable=range(len(transcripts))):
            transcript = transcripts[t_i]
            for t_j in range(t_i, len(transcripts)):
                another_transcript = transcripts[t_j]
                if transcript.overlaps(another_transcript) and transcript.gene_id != another_transcript.gene_id:
                    writer.write_feature(transcript.get_data())
                    writer.write_feature(another_transcript.get_data())
                    writer.write_comment("")


def main(args: List[str]):
    args = _parse_args(args)
    if args.out is None:
        out_basename = args.gtf
    else:
        out_basename = args.out
    gv = GeneViewFactory.from_file(args.gtf, not_save_index=True)

    with get_writer(f"{out_basename}.gene.tsv") as gene_writer, \
            get_writer(f"{out_basename}.transcripts.tsv") as transcripts_writer, \
            get_writer(f"{out_basename}.exons.tsv") as exons_writer:
        gene_writer.write("\t".join((
            "GENE_ID",
            "TRANSCRIPT_NUMBER"
        )) + "\n")
        transcripts_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "GENE_ID",
            "SPAN_LENGTH",
            "TRANSCRIBED_LENGTH",
            "EXON_NUMBER"
        )) + "\n")
        exons_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "EXON_ID",
            "TRANSCRIBED_LENGTH"
        )) + "\n")

        for gene in tqdm(desc="Iterating over genes...", iterable=gv.iter_genes()):

            gene_writer.write("\t".join((
                str(gene.gene_id),
                str(gene.number_of_transcripts)
            )) + "\n")

            transcripts = list(gene.iter_transcripts())
            for t_i in range(len(transcripts)):
                transcript = transcripts[t_i]

                transcript_transcribed_length = 0
                exons = list(transcript.iter_exons())
                for e_i in range(len(exons)):
                    exon = exons[e_i]
                    exon_length = exon.end - exon.start
                    transcript_transcribed_length += exon_length
                    exons_writer.write("\t".join((
                        transcript.transcript_id,
                        str(e_i),
                        str(exon_length)
                    )) + "\n")
                transcripts_writer.write("\t".join((
                    transcript.transcript_id,
                    transcript.gene_id,
                    str(transcript.end - transcript.start),
                    str(transcript_transcribed_length),
                    str(transcript.number_of_exons)
                )) + "\n")


    if args.get_overlapping_transcripts_from_different_gene:
        get_overlapping_transcripts_from_different_gene(out_basename, gv)

