from typing import List, Tuple, Iterable

from bioutils.datastructure.gene_view import GeneViewType
from bioutils.datastructure.gene_view_proxy import Transcript
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def assert_splice_site_existence(
        this_splice_site: List[Tuple[int, int]],
        all_splice_sites: List[List[Tuple[int, int]]]
) -> bool:
    if this_splice_site in all_splice_sites:
        return True
    else:
        all_splice_sites.append(this_splice_site)
        return False

def get_duplicated_transcript_ids(
        transcripts: Iterable[Transcript],
        by_splice_site: bool = True
) -> Iterable[str]:
    all_splice_sites: List[List[Tuple[int, int]]] = []

    if by_splice_site:
        for transcript in transcripts:
            this_splice_site = list(transcript.splice_sites)
            if assert_splice_site_existence(this_splice_site, all_splice_sites):
                lh.warn(f"Will remove {transcript.transcript_id}")
                yield transcript.transcript_id
    else:
        for transcript in transcripts:
            this_splice_site = list(transcript.exon_boundaries)
            if assert_splice_site_existence(this_splice_site, all_splice_sites):
                lh.warn(f"Will remove {transcript.transcript_id}")
                yield transcript.transcript_id



def gv_dedup(
        gv: GeneViewType,
        by_splice_site: bool = True,
        assume_no_cross_gene_duplication: bool = True
):
    """
    Remove duplicates in ``transcripts``.

    :param by_splice_site: Detect by splice sites rather than exon boundaries
    :param assume_no_cross_gene_duplication: Whether they may be duplications among genes.
    """
    lh.info("Finding transcript duplicates in gv...")
    if assume_no_cross_gene_duplication:
        transcript_ids_to_del = []
        for gene in tqdm(iterable=gv.genes.values()):
            transcript_ids_to_del.extend(get_duplicated_transcript_ids(
                transcripts=gene.transcripts.values(),
                by_splice_site=by_splice_site
            ))
    else:
        transcript_ids_to_del = list(get_duplicated_transcript_ids(
            transcripts=gv.transcripts.values(),
            by_splice_site=by_splice_site
        ))
    lh.info(f"Removing {len(transcript_ids_to_del)} transcript duplicate(s) in gv...")
    for transcript_id_to_del in transcript_ids_to_del:
        gv.del_transcript(transcript_id_to_del)
    lh.info("Removing transcript duplicate(s) in gv FIN")
