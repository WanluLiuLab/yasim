"""
as_events.py -- Generate AS Events


"""
import copy
import math
import uuid
import warnings
from typing import List, Tuple, Callable, Union
import random

from bioutils.datastructure.gene_view import GeneViewType
from bioutils.datastructure.gene_view_proxy import Transcript, Gene
from bioutils.datastructure.gv_helper import assert_splice_site_existence
from commonutils.importer.tqdm_importer import tqdm


def is_exon_skipping_able_transcript(transcript: Transcript) -> bool:
    return len(transcript.exons) >= 2


def is_intron_retention_able_transcript(transcript: Transcript) -> bool:
    return len(transcript.exons) >= 2


def generate_new_transcript_id(gene_id: str) -> str:
    return gene_id + str(uuid.uuid4())


def generate_new_transcript(transcript: Transcript) -> Transcript:
    new_transcript = copy.deepcopy(transcript)
    new_transcript.attribute['reference_transcript_id'] = new_transcript.transcript_id
    new_transcript.transcript_id = generate_new_transcript_id(transcript.gene_id)
    return new_transcript


def register_new_transcript(gv: GeneViewType, transcript: Transcript):
    gv.genes[transcript.gene_id].transcripts[transcript.transcript_id] = transcript
    gv.transcripts[transcript.transcript_id] = transcript


def perform_exon_skipping(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # set the percent of knocked out exons
    percent = 0.01*(random.randint(1,30))
    # get the number of exons n to be knocked out by multiplying total exon number of transcript
    trans_len = len(new_transcript.exons)
    exonKO = math.ceil(trans_len * percent)
    # randomly delete n exons from the transcript
    exon_keep = trans_len - exonKO
    new_transcript.exons = random.sample(new_transcript.exons, exon_keep)
    # refresh the exon list
    new_transcript.sort_exons()
    return new_transcript


def perform_intron_retention(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # randomly pick an exon for retention
    exon_num = random.randint(0,(len(new_transcript.exons)-2))
    # get the end coordinate of the neighbour exon
    end_pos = new_transcript.exons[exon_num+1].end
    # merge the coordinates of the two exons as the coordinate of the first exon
    new_transcript.exons[exon_num].end = end_pos
    # delete the neighbour exon
    del new_transcript.exons[exon_num+1]
    # refresh the exon list
    new_transcript.sort_exons()
    return new_transcript


def perform_alternative_3p_splicing(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # randomly pick an exon for splicing
    exon_num = random.randint(0, (len(new_transcript.exons) - 1))
    # randomly generate the percent to shorten
    exon_len = new_transcript.exons[exon_num].end - new_transcript.exons[exon_num].start
    splice_perc = 0.01 * (random.randint(1, 30))
    splice_len = math.ceil(exon_len * splice_perc)
    # change the end coordinate of the exon
    new_transcript.exons[exon_num].end -= splice_len
    return new_transcript


def perform_alternative_5p_splicing(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # randomly pick an exon for splicing
    exon_num = random.randint(0, (len(new_transcript.exons) - 1))
    # randomly generate the percent to shorten
    exon_len = new_transcript.exons[exon_num].end - new_transcript.exons[exon_num].start
    splice_perc = 0.01 * (random.randint(1, 30))
    splice_len = math.ceil(exon_len * splice_perc)
    # change the start coordinate of the exon
    new_transcript.exons[exon_num].start -= splice_len
    return new_transcript


class ASManipulator:
    _gv: GeneViewType
    """
    Underlying GeneView
    """

    def __init__(self, gv: GeneViewType):
        self._gv = gv

    def determine_transcript_type(self) -> Callable[[Transcript], Transcript]:
        return random.choices(
            population = (
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


    def try_generate_n_isoform_for_a_gene(self, gene:Gene, n:int):
        if len(gene.transcripts) == n:
            return
        elif len(gene.transcripts) > n:
            gek = list(gene.transcripts.keys())
            while len(gene.transcripts) > n:
                gene.transcripts.pop(gek.pop())
            return
        elif len(gene.transcripts) < n:
            all_splice_sites: List[List[Tuple[int, int]]] = []
            gtv = list(gene.transcripts.values())
            while len(gene.transcripts) < n:
                number_of_fail = 0
                new_transcript = self.determine_transcript_type()(random.choice(gtv))
                this_splice_site=list(new_transcript.splice_sites)
                if not assert_splice_site_existence(this_splice_site, all_splice_sites):
                    self._gv.add_transcript(new_transcript)
                else:
                    number_of_fail += 1
                if number_of_fail > 2 * n:
                    raise ValueError("Generation FAILED!")
                elif len(gene.transcripts) == n:
                    return

    def run(self, mu:Union[int, float]):
        gene_ids_to_del:List[str] = []
        for gene in tqdm(iterable=self._gv.genes.values(), desc="Getting exon superset..."):
            n = 0
            while n <= 0 or n >= mu*2:
                n = int(random.gauss(mu, 1))
            try:
                self.try_generate_n_isoform_for_a_gene(gene, n)
            except ValueError:
                gene_ids_to_del.append(gene.gene_id)
        for gene_id in gene_ids_to_del:
            self._gv.del_gene(gene_id)

