from typing import List, Tuple

from bioutils.datastructure.gene_view import GeneViewType
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def dedup(gv: GeneViewType, by_splice_site: bool = True):
    """
    Remove duplicates in ``transcripts``.

    :param by_splice_site: Detect by splice sites rather than exon boundaries
    """
    lh.info("Finding transcript duplicates in gv...")

    def assert_existence():
        if this_splice_site in all_splice_sites:
            transcript_to_del.append(transcript.transcript_id)
            lh.warn(f"Will remove {transcript.transcript_id}")
        else:
            all_splice_sites.append(this_splice_site)

    transcript_to_del = []
    all_splice_sites: List[List[Tuple[int, int]]] = []
    if by_splice_site:
        for transcript in tqdm(iterable=gv.transcripts.values(), desc="Iterating transcripts"):
            this_splice_site = list(transcript.splice_sites)
            assert_existence()
    else:
        for transcript in tqdm(iterable=gv.transcripts.values(), desc="Iterating transcripts"):
            this_splice_site = list(transcript.exon_boundaries)
            assert_existence()
    lh.info("Removing transcript duplicates in gv...")
    map(gv.del_transcript, transcript_to_del)
    lh.info("Removing transcript duplicates in gv FIN")
