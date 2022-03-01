"""Get satistics about GTF files that can be parsed into a Gene-Transcript-Exon Three-Tier Structure"""

import argparse
import statistics
from typing import List

from matplotlib import pyplot as plt

from bioutils.datastructure.gene_view import GeneView
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args(args)


def stat(item: List[int], fig_name: str):
    plt.clf()
    quantiles = statistics.quantiles(item)
    print(f"{fig_name} " + ", ".join((
        f"min={min(item)}",
        f"mean={round(statistics.mean(item), 2)}",
        f"median={statistics.median(item)}",
        f"max={max(item)}",
        f"quantiles={quantiles}"
    )))
    plt.hist(item, bins=150)
    # plt.xlim(0, quantiles[2])
    plt.savefig(f"{fig_name}_distribution.png")


def main(args: List[str]):
    args = _parse_args(args)
    out_basename = args.out
    gv = GeneView.from_file(args.gtf)
    # transcript_numbers = []
    # """Total number of transcripts"""
    #
    # exon_numbers = []
    # exon_length = []
    # start_and_end_sites_number = []
    # gene_with_antisense_transcripts = defaultdict(lambda: [])
    # gene_with_antisense_transcripts_on_same_chr = defaultdict(lambda: [])
    # transcript_with_antisense_exons = defaultdict(lambda: [])

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

    for gene in tqdm(desc="Iterating over genes...", iterable=gv.genes.values()):

        gene_writer.write("\t".join((
            str(gene.gene_id),
            str(len(gene.transcripts))
        )) + "\n")

        transcripts = list(gene.transcripts.values())
        for t_i in range(len(transcripts)):
            transcript = transcripts[t_i]

            transcript_transcribed_length = 0
            exons = list(transcript.exons)
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
                len(transcript.exons)
            )) + "\n")

            for t_j in range(t_i, len(transcripts)):
                another_transcript = transcripts[t_j]
                if transcript.strand != another_transcript.strand:
                    gene_with_antisense_transcripts[gene.gene_id].append(
                        f"{transcript.to_gtf_record()}\n{another_transcript.to_gtf_record()}\n\n"
                    )
                    if transcript.seqname == another_transcript.seqname:
                        gene_with_antisense_transcripts_on_same_chr[gene.gene_id].append(
                            f"{transcript.to_gtf_record()}\n{another_transcript.to_gtf_record()}\n\n"
                        )

    transcripts = list(gv.transcripts.values())
    with get_writer("overlapping_transcript.gtf") as writer:
        for t_i in tqdm(desc="Iterating over transcripts...", iterable=range(len(transcripts))):
            transcript = transcripts[t_i]
            for t_j in range(t_i, len(transcripts)):
                another_transcript = transcripts[t_j]
                if transcript.overlaps(another_transcript) and transcript.gene_id != another_transcript.gene_id:
                    writer.write(f"{transcript.to_gtf_record()}\n{another_transcript.to_gtf_record()}\n\n")

    print(f"gene_with_antisense_transcripts: {len(gene_with_antisense_transcripts)}")
    with get_writer("gene_with_antisense_transcripts.gtf") as writer:
        for gtf_str in gene_with_antisense_transcripts.values():
            writer.writelines(gtf_str)
    print(f"transcript_with_antisense_exons: {len(transcript_with_antisense_exons)}")
    with get_writer("transcript_with_antisense_exons.gtf") as writer:
        for gtf_str in transcript_with_antisense_exons.values():
            writer.writelines(gtf_str)
    print(f"gene_with_antisense_transcripts_on_same_chr: {len(gene_with_antisense_transcripts_on_same_chr)}")
    with get_writer("gene_with_antisense_transcripts_on_same_chr.gtf") as writer:
        for gtf_str in gene_with_antisense_transcripts_on_same_chr.values():
            writer.writelines(gtf_str)