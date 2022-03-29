"""
as_events.py -- Generate AS Events


"""
import copy
import uuid
from typing import List

from bioutils.datastructure.gene_view import GeneViewFactory, GeneViewType
from bioutils.datastructure.gene_view_proxy import Transcript


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
    # TODO
    return new_transcript


def perform_intron_retention(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # TODO
    return new_transcript


def perform_alternative_3p_splicing(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # TODO
    return new_transcript


def perform_alternative_5p_splicing(transcript: Transcript) -> Transcript:
    new_transcript = generate_new_transcript(transcript)
    # TODO
    return new_transcript


class ASManipulator:
    _gv: GeneViewType
    """
    Underlying GeneView
    """

    intron_retentionable_transcript: List[Transcript]
    """
    Transcript that may have intron retention events.
    """

    def __init__(self, gv: GeneViewFactory):
        self._gv = gv

    def determine_AS_able_transcripts(self):
        pass

    def simulate_as_events(self):
        pass
