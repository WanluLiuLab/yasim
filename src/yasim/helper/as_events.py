"""
as_events.py -- Generate AS Events

This file contains necessary modules in *de novo* generation of AS events.
It would reshape Number of Isoforms per Gene on Reference Genome in a log-normal distribution,
with new isoforms introduced by creating AS events and redundant isoforms removed.
"""

# FIXME: May exceed chromosome bound!

__all__ = (
    "ASManipulator",
)

import random
from typing import List, Callable, Iterable

import numpy as np
from scipy.stats import lognorm

from labw_utils.bioutils.datastructure.gene_view import GeneViewType
from labw_utils.bioutils.datastructure.gv_feature_proxy import Gene
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

ORGANISM_PARAMS = {
    "ce": (1.0274021145895147, 0.6524307003952217, 0.41902707818534024)
}
"""
(s, loc, scale) of a log-normal distribution.
"""

ORGANISM_WEIGHTS = {
    "ce": (
        11714,
        5886,
        8222,
        4371
    )
}
"""
Weight of each different AS event type.
"""

_lh = get_logger(__name__)


class ImpossibleToGenerateASEventError(ValueError):
    ...


class ASManipulator:
    _gv: GeneViewType

    def __init__(self, gv: GeneViewType):
        """
        Creation of an AS manipulator.

        :param gv: GeneView.

        .. warning::
            The Gene View passed would be internally modified!
        """
        self._gv = gv
        for transcript_id in self._gv.iter_transcript_ids():
            self._gv.transcript_sort_exons(transcript_id, "none")

    def _generate_one_as_event(
            self,
            core_func: Callable[[str], str],
            gene: Gene,
            max_try: int
    ):
        for _ in range(max_try):
            transcript_id = random.choice(list(gene.iter_transcript_ids()))
            try:
                new_transcript_id = core_func(transcript_id)
            except IndexError:
                # lh.debug(f"{core_func.__name__} IndexError")
                continue
            self._gv.transcript_sort_exons(new_transcript_id, "none")
            duplicated_transcript_id = gene.check_whether_one_transcript_duplicates_with_others(new_transcript_id)
            if duplicated_transcript_id is not None:
                self._gv.del_transcript(new_transcript_id)
                # lh.debug(f"{core_func.__name__} Duplicated")
            elif gene.get_transcript(new_transcript_id).transcribed_length < 250:
                self._gv.del_transcript(new_transcript_id)
                # lh.debug(f"{core_func.__name__} Too short")
            else:
                return
        _lh.debug(f"{gene.gene_id} trials exhausted")
        raise ImpossibleToGenerateASEventError

    def _generate_multiple_as_events(
            self,
            core_funcs: Iterable[Callable[[str], str]],
            gene: Gene,
            max_try: int = 100
    ):
        for core_func in core_funcs:
            self._generate_one_as_event(core_func=core_func, gene=gene, max_try=max_try)

    def _core_perform_exon_skipping(self, transcript_id: str) -> str:
        number_of_exons = self._gv.get_transcript(transcript_id).number_of_exons
        # set the percent of knocked out exons
        percent = random.random() * 0.3 + 0.3
        # get the number of exons n to be knocked out by multiplying total exon number of transcript
        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        es_ids = random.sample(range(0, number_of_exons), int(number_of_exons * percent))
        if len(es_ids) == 0:
            raise ImpossibleToGenerateASEventError
        _lh.debug(f"{new_transcript_id}: ES {es_ids}, total={number_of_exons}")
        for exon_id in es_ids:
            self._gv.del_exon(new_transcript_id, exon_id)
        return new_transcript_id

    def _core_perform_intron_retention(self, transcript_id: str) -> str:
        number_of_exons = self._gv.get_transcript(transcript_id).number_of_exons
        start_exon_id = random.choice(range(0, number_of_exons - 2))
        stop_exon_id = start_exon_id + 1
        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        _lh.debug(f"{new_transcript_id}: IR {start_exon_id}-{stop_exon_id}, total={number_of_exons}")
        new_transcript = self._gv.get_transcript(new_transcript_id)
        new_transcript.get_nth_exon(start_exon_id).end = new_transcript.get_nth_exon(stop_exon_id).end
        self._gv.del_exon(new_transcript_id, stop_exon_id)
        return new_transcript_id

    def _core_perform_alternative_3p_splicing(self, transcript_id: str) -> str:
        transcript = self._gv.get_transcript(transcript_id)
        # randomly pick an exon for splicing
        exon_id = random.randint(0, (transcript.number_of_exons - 1))
        # randomly generate the percent to shorten
        exon_len = transcript.get_nth_exon(exon_id).transcribed_length
        intron_len = transcript.get_intron_length(exon_id)
        splice_perc = (random.random() - 0.5) * 1.6
        delta = int(min(exon_len, intron_len) * splice_perc)
        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        _lh.debug(f"{new_transcript_id}: A3P {exon_id} (EXON_LEN={exon_len}, INTRON_LEN={intron_len}), {delta}")
        new_transcript = self._gv.get_transcript(new_transcript_id)
        new_transcript.get_nth_exon(exon_id).end += delta
        return new_transcript_id

    def _core_perform_alternative_5p_splicing(self, transcript_id: str) -> str:
        transcript = self._gv.get_transcript(transcript_id)
        # randomly pick an exon for splicing
        exon_id = random.randint(0, transcript.number_of_exons - 1)
        # randomly generate the percent to shorten
        exon_len = transcript.get_nth_exon(exon_id).transcribed_length
        intron_len = transcript.get_intron_length(exon_id - 1)
        splice_perc = (random.random() - 0.5) * 1.6
        delta = int(min(exon_len, intron_len) * splice_perc)
        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        _lh.debug(f"{new_transcript_id}: A5P {exon_id} (EXON_LEN={exon_len}, INTRON_LEN={intron_len}), {delta}")
        new_transcript = self._gv.get_transcript(new_transcript_id)
        new_transcript.get_nth_exon(exon_id).start += delta
        return new_transcript_id

    def _perform_alternative_splicing(self, gene: Gene, organism: str) -> None:
        core_funcs: Iterable[Callable[[str], str]] = random.choices(
            population=(
                (self._core_perform_alternative_3p_splicing,),
                (self._core_perform_alternative_5p_splicing,),
                (self._core_perform_exon_skipping,),
                (self._core_perform_intron_retention,),
            ),
            weights=ORGANISM_WEIGHTS[organism],
            k=1
        )[0]
        self._generate_multiple_as_events(core_funcs, gene)

    def _try_generate_n_isoform_for_a_gene(self, gene: Gene, n: int, organism: str):
        transcript_ids_to_del = []
        for transcript in gene.iter_transcripts():
            if transcript.transcribed_length < 250:
                transcript_ids_to_del.append(transcript.transcript_id)
        for transcript_id in transcript_ids_to_del:
            self._gv.del_transcript(transcript_id, auto_remove_empty_gene=False)
        if gene.number_of_transcripts == 0:
            raise ImpossibleToGenerateASEventError("Gene cannot transcribed to length >= 250")
        elif gene.number_of_transcripts == n:
            return
        elif gene.number_of_transcripts > n:
            gek = list(gene.iter_transcript_ids())
            while gene.number_of_transcripts > n:
                self._gv.del_transcript(gek.pop())
            return
        elif gene.number_of_transcripts < n:
            number_of_fail = 0
            while True:
                try:
                    self._perform_alternative_splicing(gene, organism)
                except ImpossibleToGenerateASEventError:
                    number_of_fail += 1
                if number_of_fail > 2 * n:
                    raise ImpossibleToGenerateASEventError("Generation FAILED!")
                elif gene.number_of_transcripts == n:
                    return

    def run(self, organism: str, compl_idx: int) -> None:
        """
        Execute the ASManipulator.

        :param organism: Name of the organism. Allowed are:

            - ``ce``: *C. Elegans*.
        :param compl_idx: The complexity index. Larger for larger Number of Isoforms in a Gene.
        """
        if organism not in ORGANISM_PARAMS.keys():
            raise ValueError(f"no organism {organism} available! Available: {ORGANISM_PARAMS.keys()}")
        fit_tuple = ORGANISM_PARAMS[organism]
        gene_ids_to_del: List[str] = []
        all_gene_ids = list(self._gv.iter_gene_ids())
        targeted_nipg = np.array(
            lognorm.rvs(
                fit_tuple[0],
                loc=fit_tuple[1],
                scale=fit_tuple[2],
                size=len(all_gene_ids) * 5
            ),
            dtype=float
        )
        minvalue = np.min(targeted_nipg)
        targeted_nipg = (targeted_nipg - minvalue) * compl_idx / np.mean(targeted_nipg) + minvalue
        targeted_nipg = targeted_nipg[np.logical_and(1 < targeted_nipg, targeted_nipg < 25)][:len(all_gene_ids)]
        targeted_nipg = np.array(np.sort(targeted_nipg), dtype=int)
        reference_nipg = [
            self._gv.get_gene(gene_id).number_of_transcripts for gene_id in all_gene_ids
        ]
        reference_nipg_rank = np.argsort(reference_nipg)
        for gene_index in tqdm(
                iterable=range(len(all_gene_ids)),
                desc="Generating isoforms...",
                total=self._gv.number_of_genes
        ):
            gene = self._gv.get_gene(all_gene_ids[gene_index])
            n = targeted_nipg[reference_nipg_rank[gene_index]]
            try:
                self._try_generate_n_isoform_for_a_gene(gene, n, organism)
            except ImpossibleToGenerateASEventError:
                gene_ids_to_del.append(gene.gene_id)

        _lh.info(f"Will remove {len(gene_ids_to_del)} genes out of {self._gv.number_of_genes}")
        for gene_id in gene_ids_to_del:
            self._gv.del_gene(gene_id)
        _lh.info("Will remove genes FIN")
        self._gv.standardize()

    def to_file(self, dst_gtf_file_path: str) -> None:
        """
        Write the internal GeneView to GTF.

        :param dst_gtf_file_path: Path to GTF to be written.
        """
        self._gv.to_file(dst_gtf_file_path)
