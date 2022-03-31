"""
as_events.py -- Generate AS Events


"""
import copy
import math
import uuid
from typing import List
import random

from bioutils.datastructure.gene_view import GeneViewType
from bioutils.datastructure.gene_view_proxy import Transcript
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

    def generate_exon_superset(self):
        for gene in tqdm(iterable=self._gv.genes.values(), desc="Getting exon superset..."):
            gene.generate_exon_superset()

    def run(self):
        pass
