"""
as_events.py -- Generate AS Events

This file contains necessary modules in *de novo* generation of AS events.
It would reshape Number of Isoforms per Gene on Reference Genome in a log-normal distribution,
with new isoforms introduced by creating AS events and redundant isoforms removed.

.. versionadded:: 3.1.5
"""
from __future__ import annotations

__all__ = ("ASManipulator",)

import array
import random
import uuid

import numpy as np
from scipy.stats import lognorm

from labw_utils.bioutils.datastructure.gene_tree import GeneTreeInterface
from labw_utils.bioutils.datastructure.gv.feature_proxy import BaseFeatureProxy
from labw_utils.bioutils.datastructure.gv.gene import Gene
from labw_utils.bioutils.datastructure.gv.transcript import Transcript
from labw_utils.bioutils.parser.gtf import GtfIteratorWriter
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Callable, Iterable, Set

ORGANISM_PARAMS = {"ce": (1.0274021145895147, 0.6524307003952217, 0.41902707818534024)}
"""
(s, loc, scale) of a log-normal distribution.

.. versionadded:: 3.1.5
"""

ORGANISM_WEIGHTS = {"ce": (11714, 5886, 8222, 4371)}
"""
Weight of each different AS event type.

.. versionadded:: 3.1.5
"""

_lh = get_logger(__name__)


class SpaceOptimizedTranscript(Transcript):
    __slots__ = ["_exons_bytes"]
    _exons_bytes: bytes

    def __init__(self, t: Transcript):
        self._cdna = None
        self._cdna_unspliced = None
        self._exon_boundaries = None
        self._cdna_unspliced_masked = None
        self._splice_sites = None
        self._is_inferred = t.is_inferred
        self._transcript_id = t.transcript_id  # type: ignore
        self._gene_id = t.gene_id  # type: ignore
        self._exons = list(t.exons)  # Need to copy. # type: ignore
        BaseFeatureProxy.__init__(self, data=t.get_data(), is_checked=t.is_checked)
        self.rescale_from_exon_boundaries(force=True)
        selist = []
        for exon in self.exons:
            selist.append(exon.start0b)
            selist.append(exon.end0b)
        self._exons_bytes = b"".join(
            (
                array.array("Q", selist).tobytes(),
                bytes(self.seqname, encoding="UTF-8"),
                bytes(self.strand),
            )
        )

    def __hash__(self):
        return hash(self._exon_boundaries)

    def __eq__(self, other: "SpaceOptimizedTranscript") -> bool:
        return self._exons_bytes == other._exons_bytes


class ImpossibleToGenerateASEventError(ValueError):
    """
    TODO docs

    .. versionadded:: 3.1.5
    """

    ...


class ASManipulator:
    """
    TODO docs

    .. versionadded:: 3.1.5
    """

    _gv: GeneTreeInterface
    _n_new_isoforms: int
    _n_removed_isoforms: int

    def __init__(self, gv: GeneTreeInterface):
        """
        Creation of an AS manipulator.

        :param gv: GeneTree.

        .. warning::
            The Gene View passed would be internally modified!
        """
        self._gv = gv
        self._n_new_isoforms = 0
        self._n_removed_isoforms = 0

    # @profile
    def _generate_one_as_event(
        self,
        core_func: Callable[[Transcript], Transcript],
        gene: Gene,
        max_try: int,
        sotes: Set[SpaceOptimizedTranscript],
    ) -> Gene:
        for _ in range(max_try):
            transcript = random.choice(list(gene.transcript_values))
            try:
                new_transcript = core_func(transcript)
            except IndexError:
                # lh.debug(f"{core_func.__name__} IndexError")
                continue
            if new_transcript.transcribed_length < 250:
                continue
            sotet = SpaceOptimizedTranscript(new_transcript)
            if sotet in sotes:
                continue
            gene = gene.add_transcript(new_transcript)
            sotes.add(sotet)
            self._n_new_isoforms += 1
            return gene
        _lh.debug("GEN AS EVENT: Gene %s: %d trials exhausted", gene.gene_id, max_try)
        raise ImpossibleToGenerateASEventError

    def _generate_multiple_as_events(
        self,
        core_funcs: Iterable[Callable[[Transcript], Transcript]],
        gene: Gene,
        max_try: int = 100,
    ) -> Gene:
        sotes = set(
            map(
                lambda t: SpaceOptimizedTranscript(t),
                gene.transcript_values,
            )
        )
        for core_func in core_funcs:
            gene = self._generate_one_as_event(core_func=core_func, gene=gene, max_try=max_try, sotes=sotes)
        return gene

    def _core_perform_exon_skipping(self, old_transcript: Transcript) -> Transcript:
        number_of_exons = old_transcript.number_of_exons
        # set the percent of knocked out exons
        percent = random.random() * 0.3 + 0.3
        # get the number of exons n to be knocked out by multiplying total exon number of transcript
        new_transcript = old_transcript.duplicate().update_attribute(
            transcript_id=old_transcript.gene_id + str(uuid.uuid4())
        )
        es_ids = random.sample(range(0, number_of_exons), int(number_of_exons * percent))
        if len(es_ids) == 0:
            raise ImpossibleToGenerateASEventError
        _lh.debug(f"GEN AS EVENT: {new_transcript.transcript_id}: ES {es_ids}, total={number_of_exons}")
        for exon_id in es_ids:
            new_transcript = new_transcript.del_exon(exon_id)
        return new_transcript

    def _core_perform_intron_retention(self, old_transcript: Transcript) -> Transcript:
        number_of_exons = old_transcript.number_of_exons
        start_exon_id = random.choice(range(0, number_of_exons - 2))
        stop_exon_id = start_exon_id + 1
        new_transcript = old_transcript.duplicate().update_attribute(
            transcript_id=old_transcript.gene_id + str(uuid.uuid4())
        )
        _lh.debug(
            f"GEN AS EVENT: {new_transcript.transcript_id}: IR {start_exon_id}-{stop_exon_id}, total={number_of_exons}"
        )
        exon_to_change = new_transcript.get_exon(start_exon_id)
        new_transcript = new_transcript.del_exon(start_exon_id)
        exon_to_change = exon_to_change.update_attribute(end=new_transcript.get_exon(stop_exon_id).end)
        new_transcript.add_exon(exon_to_change)
        new_transcript = new_transcript.del_exon(stop_exon_id)
        return new_transcript

    def _core_perform_alternative_3p_splicing(self, old_transcript: Transcript) -> Transcript:
        # randomly pick an exon for splicing
        exon_id = random.randint(0, (old_transcript.number_of_exons - 1))
        # randomly generate the percent to shorten
        exon_len = old_transcript.get_exon(exon_id).transcribed_length
        intron_len = old_transcript.get_intron_length(exon_id)
        splice_perc = (random.random() - 0.5) * 1.6
        delta = int(min(exon_len, intron_len) * splice_perc)
        new_transcript = old_transcript.duplicate().update_attribute(
            transcript_id=old_transcript.gene_id + str(uuid.uuid4())
        )
        _lh.debug(
            f"GEN AS EVENT: {new_transcript.transcript_id}: "
            f"A3P {exon_id} (EXON_LEN={exon_len}, INTRON_LEN={intron_len}), {delta}"
        )
        exon_to_change = new_transcript.get_exon(exon_id)
        new_transcript = new_transcript.del_exon(exon_id)
        exon_to_change = exon_to_change.update_attribute(end=exon_to_change.end + delta)
        new_transcript.add_exon(exon_to_change)
        return new_transcript

    def _core_perform_alternative_5p_splicing(self, old_transcript: Transcript) -> Transcript:
        # randomly pick an exon for splicing
        exon_id = random.randint(0, old_transcript.number_of_exons - 1)
        # randomly generate the percent to shorten
        exon_len = old_transcript.get_exon(exon_id).transcribed_length
        intron_len = old_transcript.get_intron_length(exon_id - 1)
        splice_perc = (random.random() - 0.5) * 1.6
        delta = int(min(exon_len, intron_len) * splice_perc)
        new_transcript = old_transcript.duplicate().update_attribute(
            transcript_id=old_transcript.gene_id + str(uuid.uuid4())
        )
        _lh.debug(
            f"GEN AS EVENT: {new_transcript.transcript_id}: "
            f"A5P {exon_id} (EXON_LEN={exon_len}, INTRON_LEN={intron_len}), {delta}"
        )
        exon_to_change = new_transcript.get_exon(exon_id)
        new_transcript = new_transcript.del_exon(exon_id)
        exon_to_change = exon_to_change.update_attribute(start=exon_to_change.start + delta)
        new_transcript.add_exon(exon_to_change)
        return new_transcript

    def _perform_alternative_splicing(self, gene: Gene, organism: str) -> Gene:
        core_funcs: Iterable[Callable[[str], SpaceOptimizedTranscript]] = random.choices(
            population=(
                (self._core_perform_alternative_3p_splicing,),
                (self._core_perform_alternative_5p_splicing,),
                (self._core_perform_exon_skipping,),
                (self._core_perform_intron_retention,),
            ),
            weights=ORGANISM_WEIGHTS[organism],
            k=1,
        )[0]
        return self._generate_multiple_as_events(core_funcs, gene)

    def _try_generate_n_isoform_for_a_gene(self, gene: Gene, n: int, organism: str):
        transcript_ids_to_del = []
        for transcript in gene.transcript_values:
            if transcript.transcribed_length < 250:
                transcript_ids_to_del.append(transcript.transcript_id)
        for transcript_id in transcript_ids_to_del:
            gene = gene.del_transcript(transcript_id)
            self._n_removed_isoforms += 1
        if gene.number_of_transcripts == 0:
            raise ImpossibleToGenerateASEventError("Gene cannot transcribed to length >= 250")
        elif gene.number_of_transcripts == n:
            pass
        elif gene.number_of_transcripts > n:
            while gene.number_of_transcripts > n:
                gene = gene.del_transcript(gene.transcript_ids[0])
                self._n_removed_isoforms += 1
        elif gene.number_of_transcripts < n:
            number_of_fail = 0
            while True:
                try:
                    gene = self._perform_alternative_splicing(gene, organism)
                except ImpossibleToGenerateASEventError:
                    number_of_fail += 1
                if number_of_fail > 2 * n:
                    raise ImpossibleToGenerateASEventError("Generation FAILED!")
                elif gene.number_of_transcripts == n:
                    break
        self._gv = self._gv.replace_gene(gene)

    def run(self, organism: str, compl_idx: int) -> None:
        """
        Execute the ASManipulator.

        :param organism: Name of the organism. Allowed are:

            - ``ce``: *C. Elegans*.
        :param compl_idx: The Transcriptome Complexity Index. Larger for larger Number of Isoforms in a Gene.
        """

        _lh.info(
            "GEN AS EVENT: Loaded %d genes with %d transcript",
            self._gv.number_of_genes,
            self._gv.number_of_transcripts,
        )
        if organism not in ORGANISM_PARAMS.keys():
            raise ValueError(f"no organism {organism} available! Available: {ORGANISM_PARAMS.keys()}")
        fit_tuple = ORGANISM_PARAMS[organism]
        gene_ids_to_del: List[str] = []
        all_gene_ids = list(self._gv.gene_ids)
        targeted_nipg = np.array(
            lognorm.rvs(fit_tuple[0], loc=fit_tuple[1], scale=fit_tuple[2], size=len(all_gene_ids) * 5), dtype=float
        )
        minvalue = np.min(targeted_nipg)
        targeted_nipg = (targeted_nipg - minvalue) * compl_idx / np.mean(targeted_nipg) + minvalue
        targeted_nipg = targeted_nipg[np.logical_and(1 < targeted_nipg, targeted_nipg < 25)][: len(all_gene_ids)]
        targeted_nipg = np.array(np.sort(targeted_nipg), dtype=int)
        reference_nipg = [self._gv.get_gene(gene_id)[0].number_of_transcripts for gene_id in all_gene_ids]
        reference_nipg_rank = np.argsort(reference_nipg)
        for gene_index in tqdm(
            iterable=range(len(all_gene_ids)),
            desc="Generating isoforms...",
            total=self._gv.number_of_genes,
        ):
            gene = self._gv.get_gene(all_gene_ids[gene_index])[0]
            n = targeted_nipg[reference_nipg_rank[gene_index]].item()
            try:
                self._try_generate_n_isoform_for_a_gene(gene, n, organism)
            except ImpossibleToGenerateASEventError:
                gene_ids_to_del.append(gene.gene_id)
                self._n_removed_isoforms += gene.number_of_transcripts

        _lh.info("GEN AS EVENT: Removing genes that failed to produce enough isoforms...")
        for gene_id in tqdm(gene_ids_to_del, desc="Removing genes"):
            self._gv = self._gv.del_gene(gene_id)
        _lh.info(
            "GEN AS EVENT: Finished with %d genes (-%d) and %d transcript (+%d, -%d)",
            self._gv.number_of_genes,
            len(gene_ids_to_del),
            self._gv.number_of_transcripts,
            self._n_new_isoforms,
            self._n_removed_isoforms,
        )

    def to_file(self, dst_gtf_file_path: str) -> None:
        """
        Write the internal GeneView to GTF.

        :param dst_gtf_file_path: Path to GTF to be written.
        """
        GtfIteratorWriter.write_iterator(self._gv.to_feature_iterator(), dst_gtf_file_path)
