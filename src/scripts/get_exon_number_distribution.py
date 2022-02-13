import argparse
import statistics
import sys
from collections import defaultdict
from typing import List

from bioutils.datastructure import GeneView
from matplotlib import pyplot as plt

from commonutils import ioctl
from commonutils.tqdm_importer import tqdm


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
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
    gv = GeneView.from_file(args.gtf)
    transcript_numbers = []
    gene_span_length = []
    transcript_span_length = []
    transcript_length = []
    exon_numbers = []
    exon_length = []
    start_and_end_sites_number = []
    gene_with_antisense_transcripts = defaultdict(lambda: [])
    gene_with_antisense_transcripts_on_same_chr = defaultdict(lambda: [])
    transcript_with_antisense_exons = defaultdict(lambda: [])
    overlapping_transcript = defaultdict(lambda: [])

    for gene in tqdm(desc="Iterating over genes...", iterable=gv.genes.values()):
        max_transcript_span_length = 0
        start_sites = set()
        end_sites = set()
        transcript_numbers.append(len(gene.transcripts))
        transcripts = list(gene.transcripts.values())
        for t_i in range(len(transcripts)):
            transcript = transcripts[t_i]
            exon_numbers.append(len(transcript.exons))
            tmp_transcript_length = transcript.end - transcript.start
            transcript_span_length.append(tmp_transcript_length)
            max_transcript_span_length = max(max_transcript_span_length, tmp_transcript_length)
            exons = list(transcript.exons)
            for e_i in range(len(exons)):
                exon = exons[e_i]
                tmp_exon_length = exon.end - exon.start
                tmp_transcript_length += tmp_exon_length
                exon_length.append(tmp_exon_length)
                start_sites.add((exon.seqname, exon.start))
                end_sites.add((exon.seqname, exon.end))
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
            transcript_length.append(tmp_transcript_length)
        start_and_end_sites = start_sites
        start_and_end_sites.update(end_sites)
        start_and_end_sites_number.append(len(start_and_end_sites))
        gene_span_length.append(max_transcript_span_length)

    transcripts = list(gv.transcripts.values())
    with ioctl.get_writer("overlapping_transcript.gtf") as writer:
        for t_i in tqdm(desc="Iterating over transcripts...",iterable=range(len(transcripts))):
            transcript = transcripts[t_i]
            for t_j in range(t_i, len(transcripts)):
                another_transcript = transcripts[t_j]
                if transcript.overlaps(another_transcript) and transcript.gene_id != another_transcript.gene_id:
                    writer.write(f"{transcript.to_gtf_record()}\n{another_transcript.to_gtf_record()}\n\n")

    gene_with_antisense_transcripts = gene_with_antisense_transcripts
    transcript_with_antisense_exons = transcript_with_antisense_exons
    stat(transcript_numbers, "transcript_numbers_in_a_gene")
    stat(exon_numbers, "exon_numbers_in_a_transcript")
    stat(gene_span_length, "gene_span_length")
    stat(transcript_span_length, "transcript_span_length")
    stat(exon_length, "exon_length")
    stat(transcript_length, "transcript_length")
    stat(start_and_end_sites_number, "start_and_end_sites_number_in_a_gene")
    print(f"gene_with_antisense_transcripts: {len(gene_with_antisense_transcripts)}")
    with ioctl.get_writer("gene_with_antisense_transcripts.gtf") as writer:
        for gtf_str in gene_with_antisense_transcripts.values():
            writer.writelines(gtf_str)
    print(f"transcript_with_antisense_exons: {len(transcript_with_antisense_exons)}")
    with ioctl.get_writer("transcript_with_antisense_exons.gtf") as writer:
        for gtf_str in transcript_with_antisense_exons.values():
            writer.writelines(gtf_str)
    print(f"gene_with_antisense_transcripts_on_same_chr: {len(gene_with_antisense_transcripts_on_same_chr)}")
    with ioctl.get_writer("gene_with_antisense_transcripts_on_same_chr.gtf") as writer:
        for gtf_str in gene_with_antisense_transcripts_on_same_chr.values():
            writer.writelines(gtf_str)

if __name__ == "__main__":
    main(sys.argv[1:])
