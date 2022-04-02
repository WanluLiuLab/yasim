import copy
import uuid
from typing import List, Tuple, Iterable

from bioutils.datastructure.gene_view import GeneViewType
from bioutils.datastructure.gene_view_proxy import Transcript, Gene, Exon
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


def generate_new_transcript_id(gene_id: str) -> str:
    return gene_id + str(uuid.uuid4())


def generate_new_transcript(transcript: Transcript) -> Transcript:
    new_transcript = copy.deepcopy(transcript)
    new_transcript.attribute['reference_transcript_id'] = new_transcript.transcript_id
    new_transcript.transcript_id = generate_new_transcript_id(transcript.gene_id)
    return new_transcript


def enable_exon_superset():
    Gene.register_new_attribute("_exon_superset")

    def generate_exon_superset(self: Gene):
        def add_exon(_all_exons: List[Exon], new_exon: Exon):
            for _exon in _all_exons:
                if new_exon == _exon:
                    return
            _all_exons.append(new_exon)

        if self._exon_superset is not None:
            return
        self._exon_superset: List[Exon] = []
        for transcript in self.transcripts.values():
            for exon in transcript.exons:
                add_exon(self._exon_superset, exon)

    def get_exon_superset(self: Gene):
        return self._exon_superset

    Gene.register_new_hook(generate_exon_superset)
    Gene.register_new_hook(get_exon_superset)
