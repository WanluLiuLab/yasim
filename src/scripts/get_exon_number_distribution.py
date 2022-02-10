import argparse
import statistics
import sys
from typing import List

from matplotlib import pyplot as plt

from bioutils.datastructure.gene import GeneView
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
    gv = GeneView(args.gtf)
    transcript_numbers = []
    gene_span_length = []
    transcript_span_length = []
    transcript_length = []
    exon_numbers = []
    exon_length = []
    start_and_end_sites_number = []
    gene_with_start_end_equal = []
    gene_with_antisense_transcripts = []
    transcript_with_antisense_exons = []

    for gene in tqdm(desc="Iterating over genes...", iterable=gv.genes.values()):
        max_transcript_span_length = 0
        start_sites = set()
        end_sites = set()
        transcript_numbers.append(len(gene.transcripts))
        for transcript in gene.transcripts.values():
            exon_numbers.append(len(transcript.exons))
            tmp_transcript_length = transcript.end - transcript.start
            transcript_span_length.append(tmp_transcript_length)
            max_transcript_span_length = max(max_transcript_span_length, tmp_transcript_length)
            for exon in transcript.exons:
                tmp_exon_length = exon.end - exon.start
                tmp_transcript_length += tmp_exon_length
                exon_length.append(tmp_exon_length)
                start_sites.add((exon.seqname, exon.start))
                end_sites.add((exon.seqname, exon.end))
                if exon.strand != transcript.strand:
                    transcript_with_antisense_exons.append(gene.name)
            if transcript.strand != gene.strand:
                gene_with_antisense_transcripts.append(gene.name)
            transcript_length.append(tmp_transcript_length)
        for item in start_sites:
            if item in end_sites:
                gene_with_start_end_equal.append(gene.name)
                break
        start_and_end_sites = start_sites
        start_and_end_sites.update(end_sites)
        start_and_end_sites_number.append(len(start_and_end_sites))
        gene_span_length.append(max_transcript_span_length)

    gene_with_antisense_transcripts=list(set(gene_with_antisense_transcripts))
    transcript_with_antisense_exons = list(set(transcript_with_antisense_exons))
    stat(transcript_numbers, "transcript_numbers_in_a_gene")
    stat(exon_numbers, "exon_numbers_in_a_transcript")
    stat(gene_span_length, "gene_span_length")
    stat(transcript_span_length, "transcript_span_length")
    stat(exon_length, "exon_length")
    stat(transcript_length, "transcript_length")
    stat(start_and_end_sites_number, "start_and_end_sites_number_in_a_gene")
    print(f"gene_with_start_end_equal: {len(gene_with_start_end_equal)}={gene_with_start_end_equal}")
    print(f"gene_with_antisense_transcripts: {len(gene_with_antisense_transcripts)}={gene_with_antisense_transcripts}")
    print(f"transcript_with_antisense_exons: {len(transcript_with_antisense_exons)}={transcript_with_antisense_exons}")
    pass


if __name__ == "__main__":
    main(sys.argv[1:])
