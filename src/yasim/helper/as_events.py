"""
as_events.py -- Generate AS Events


"""
import copy
import uuid
from typing import List

from bioutils.datastructure.gene_view import GeneView
from bioutils.datastructure.gene_view_proxy import Transcript


def is_exon_skipping_able_transcript(transcript: Transcript) -> bool:
    return len(transcript.exons) >= 2


def is_intron_retention_able_transcript(transcript: Transcript) -> bool:
    return len(transcript.exons) >= 2


def generate_new_transcript_id() -> str:
    return str(uuid.uuid4())


def generate_new_transcript(transcript: Transcript) -> Transcript:
    new_transcript = copy.deepcopy(transcript)
    new_transcript.attribute['reference_transcript_id'] = new_transcript.transcript_id
    new_transcript.transcript_id = generate_new_transcript_id()
    return new_transcript


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
    _gv: GeneView
    """
    Underlying GeneView
    """

    intron_retentionable_transcript: List[Transcript]
    """
    Transcript that may have intron retention events.
    """

    def __init__(self, gv: GeneView):
        self._gv = gv

    def determine_AS_able_transcripts(self):
        pass

    def simulate_as_events(self):
        pass



