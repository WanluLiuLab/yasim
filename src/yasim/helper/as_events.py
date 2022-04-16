"""
as_events.py -- Generate AS Events
"""
import math
import random
from typing import List, Tuple, Callable, Union, Iterable

from bioutils.datastructure.gene_view import GeneViewType, GeneViewFactory
from bioutils.datastructure.gene_view_proxy import Transcript, Gene
from bioutils.datastructure.gv_helper import assert_splice_site_existence, generate_new_transcript
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


# def perform_exon_skipping(
#         transcript: Transcript,
#         exon_number_to_ko: Iterable[int]
# ) -> Transcript:
#     new_transcript = generate_new_transcript(transcript)
#     sorted(exon_number_to_ko)
#
#     # set the percent of knocked out exons
#     percent = 0.01 * (random.randint(1, 30))
#     # get the number of exons n to be knocked out by multiplying total exon number of transcript
#     trans_len = len(new_transcript.exons)
#     exonKO = math.ceil(trans_len * percent)
#     # randomly delete n exons from the transcript
#     exon_keep = trans_len - exonKO
#     new_transcript.exons = random.sample(new_transcript.exons, exon_keep)
#     # refresh the exon list
#     new_transcript.sort_exons()
#     return new_transcript


def perform_exon_skipping(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # set the percent of knocked out exons
    percent = 0.01 * (random.randint(1, 30))
    # get the number of exons n to be knocked out by multiplying total exon number of transcript
    trans_len = new_transcript.number_of_exons
    exonKO = math.ceil(trans_len * percent)
    # randomly delete n exons from the transcript
    exon_keep = trans_len - exonKO
    new_transcript._exons = random.sample(new_transcript._exons, exon_keep)
    # refresh the exon list
    new_transcript.sort_exons()
    return new_transcript


def perform_intron_retention(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # randomly pick an exon for retention
    exon_num = random.randint(0, (new_transcript.number_of_exons - 2))
    # get the end coordinate of the neighbour exon
    end_pos = new_transcript._exons[exon_num + 1].end
    # merge the coordinates of the two exons as the coordinate of the first exon
    new_transcript._exons[exon_num].end = end_pos
    # delete the neighbour exon
    del new_transcript._exons[exon_num + 1]
    # refresh the exon list
    new_transcript.sort_exons()
    return new_transcript


def perform_alternative_3p_splicing(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # randomly pick an exon for splicing
    exon_num = random.randint(0, (new_transcript.number_of_exons - 1))
    # randomly generate the percent to shorten
    exon_len = new_transcript._exons[exon_num].end - new_transcript._exons[exon_num].start
    splice_perc = 0.01 * (random.randint(1, 30))
    splice_len = math.ceil(exon_len * splice_perc)
    # change the end coordinate of the exon
    new_transcript._exons[exon_num].end -= splice_len
    return new_transcript


def perform_alternative_5p_splicing(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # randomly pick an exon for splicing
    exon_num = random.randint(0, (new_transcript.number_of_exons - 1))
    # randomly generate the percent to shorten
    exon_len = new_transcript._exons[exon_num].end - new_transcript._exons[exon_num].start
    splice_perc = 0.01 * (random.randint(1, 30))
    splice_len = math.ceil(exon_len * splice_perc)
    # change the start coordinate of the exon
    new_transcript._exons[exon_num].start -= splice_len
    return new_transcript


class ASManipulator:
    _gv: GeneViewType
    """
    Underlying GeneView
    """

    def __init__(self, gv: GeneViewType):
        self._gv = gv

    def perform_alternative_splicing(self) -> Callable[[Transcript], Transcript]:
        return random.choices(
            population=(
                perform_alternative_3p_splicing,
                perform_intron_retention,
                perform_alternative_5p_splicing,
                perform_exon_skipping,
                lambda transcript: perform_alternative_3p_splicing(perform_alternative_5p_splicing(transcript)),
                lambda transcript: perform_alternative_3p_splicing(perform_exon_skipping(transcript)),
                lambda transcript: perform_alternative_3p_splicing(perform_intron_retention(transcript))
            ),
            weights=(
                832,
                717,
                568,
                333,
                173,
                118,
                111
            ),
            k=1
        )[0]

    def try_generate_n_isoform_for_a_gene(self, gene: Gene, n: int):
        transcript_ids_to_del = []
        for transcript in gene.iter_transcripts():
            if transcript.transcribed_length < 250:
                transcript_ids_to_del.append(transcript.transcript_id)
        for transcript_id in transcript_ids_to_del:
            self._gv.del_transcript(transcript_id, auto_remove_empty_gene=False)
        if gene.number_of_transcripts == 0:
            raise ValueError("Generation FAILED!")
        elif gene.number_of_transcripts == n:
            return
        elif gene.number_of_transcripts > n:
            gek = list(gene.iter_transcript_ids())
            while gene.number_of_transcripts > n:
                gene.del_transcript(gek.pop())
            return
        elif gene.number_of_transcripts < n:
            all_splice_sites: List[List[Tuple[int, int]]] = []
            while gene.number_of_transcripts < n:
                number_of_fail = 0
                new_transcript = self.perform_alternative_splicing()(random.choice(list(gene.iter_transcripts())))
                this_splice_site = list(new_transcript.splice_sites)
                if new_transcript.transcribed_length >= 250 and not assert_splice_site_existence(this_splice_site,
                                                                                                 all_splice_sites):
                    self._gv.add_transcript(new_transcript)
                else:
                    number_of_fail += 1
                if number_of_fail > 2 * n:
                    raise ValueError("Generation FAILED!")
                elif gene.number_of_transcripts == n:
                    return

    def run(self, mu: Union[int, float]):
        gene_ids_to_del: List[str] = []
        for gene in tqdm(iterable=self._gv.iter_genes(), desc="Generating isoforms..."):
            n = 0
            while n <= 0 or n >= mu * 2:
                n = int(random.gauss(mu, mu / 2))
            try:
                self.try_generate_n_isoform_for_a_gene(gene, n)
            except ValueError:
                gene_ids_to_del.append(gene.gene_id)
        lh.info(f"Will remove {len(gene_ids_to_del)} genes")
        for gene_id in gene_ids_to_del:
            self._gv.del_gene(gene_id)
        lh.info(f"Will remove genes FIN")
        self._gv.standardize()

    def to_file(self, filename: str):
        self._gv.to_file(filename)


if __name__ == '__main__':
    for i in [3]:
        gv = GeneViewFactory.from_file("/media/yuzj/BUP/iter_3/ce11.ncbiRefSeq.gtf")
        asm = ASManipulator(gv)
        asm.run(i)
        asm.to_file(f"{i}.gtf")
