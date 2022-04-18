"""
as_events.py -- Generate AS Events
"""

from __future__ import annotations

import math
import random
from typing import List, Callable, Union

from bioutils.datastructure.gene_view import GeneViewType, GeneViewFactory
from bioutils.datastructure.gene_view_proxy import Gene
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


class ImpossibleToGenerateASEventError(ValueError):
    pass


class ASManipulator:
    _gv: GeneViewType
    """
    Underlying GeneView
    """

    def __init__(self, gv: GeneViewType):
        self._gv = gv

    def core_perform_exon_skipping(self, transcript_id: str) -> str:
        number_of_exons = self._gv.get_transcript(transcript_id).number_of_exons
        # set the percent of knocked out exons
        percent = 0.01 * (random.randint(1, 30))
        # get the number of exons n to be knocked out by multiplying total exon number of transcript
        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        es_ids = random.sample(range(0, number_of_exons), math.ceil(number_of_exons * percent))
        lh.debug(f"{new_transcript_id}: ES {es_ids}")
        for exon_id in es_ids:
            self._gv.del_exon(new_transcript_id, exon_id)
        return new_transcript_id

    def general_perform_wrapper(
            self, gene: Gene,
            core_func: Callable[[str], str],
            max_try: int = 100
    ):
        transcript_ids = list(gene.iter_transcript_ids())
        transcript_id = random.choice(transcript_ids)
        for _ in range(max_try):
            try:
                new_transcript_id = core_func(transcript_id)
            except IndexError:
                lh.debug(core_func.__name__, "IndexError")
                continue
            duplicated_transcript_id = gene.check_whether_one_transcript_duplicates_with_others(new_transcript_id)
            if duplicated_transcript_id is not None:
                # print(gv.get_transcript(duplicated_transcript_id).transcript_id)
                # print(gv.get_transcript(new_transcript_id)._data.attribute['reference_transcript_id'])
                # print(gv.get_transcript(new_transcript_id).transcript_id)
                self._gv.del_transcript(new_transcript_id)
                lh.debug(core_func.__name__, "Duplicated")
            elif gene.get_transcript(new_transcript_id).transcribed_length < 250:
                self._gv.del_transcript(new_transcript_id)
                lh.debug(core_func.__name__, "Too short")
            else:
                return
        raise ImpossibleToGenerateASEventError

    def perform_exon_skipping(self, gene: Gene):
        self.general_perform_wrapper(gene, self.core_perform_exon_skipping)

    def core_perform_intron_retention(self, transcript_id: str) -> str:
        number_of_exons = self._gv.get_transcript(transcript_id).number_of_exons
        start_exon_id = random.choice(range(0, number_of_exons - 2))
        stop_exon_id = start_exon_id + 1
        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        lh.debug(f"{new_transcript_id}: IR {start_exon_id}-{stop_exon_id}")
        new_transcript = self._gv.get_transcript(new_transcript_id)
        new_transcript.get_nth_exon(start_exon_id).end = new_transcript.get_nth_exon(stop_exon_id).end
        for exon_id_to_del in range(start_exon_id + 1, stop_exon_id):
            gv.del_exon(new_transcript_id, exon_id_to_del)
        self._gv.transcript_sort_exons(new_transcript_id)
        return new_transcript_id

    def perform_intron_retention(self, gene: Gene):
        self.general_perform_wrapper(gene, self.core_perform_intron_retention)

    def core_perform_alternative_3p_splicing(self, transcript_id: str) -> str:
        transcript = self._gv.get_transcript(transcript_id)
        # randomly pick an exon for splicing
        exon_id = random.randint(0, (transcript.number_of_exons - 1))
        # randomly generate the percent to shorten
        exon_len = transcript.get_nth_exon(exon_id).end - transcript.get_nth_exon(exon_id).start
        if exon_id + 1 != transcript.number_of_exons:
            intron_len = transcript.get_nth_exon(exon_id + 1).start - transcript.get_nth_exon(exon_id).end
        else:
            intron_len = math.inf
        splice_perc = 0.01 * (random.randint(-90, 90))
        delta = math.ceil(min(exon_len, intron_len) * splice_perc)

        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        lh.debug(f"{new_transcript_id}: A3P {exon_id}, {delta}")
        new_transcript = self._gv.get_transcript(new_transcript_id)
        new_transcript.get_nth_exon(exon_id).end += delta
        self._gv.transcript_sort_exons(new_transcript_id)
        return new_transcript_id

    def perform_alternative_3p_splicing(self, gene: Gene):
        self.general_perform_wrapper(gene, self.core_perform_alternative_3p_splicing)

    def core_perform_alternative_5p_splicing(self, transcript_id: str) -> str:
        transcript = self._gv.get_transcript(transcript_id)
        # randomly pick an exon for splicing
        exon_id = random.randint(0, (transcript.number_of_exons - 1))
        # randomly generate the percent to shorten
        exon_len = transcript.get_nth_exon(exon_id).end - transcript.get_nth_exon(exon_id).start
        if exon_id != 0:
            intron_len = transcript.get_nth_exon(exon_id - 1).end - transcript.get_nth_exon(exon_id).start
        else:
            intron_len = math.inf
        splice_perc = 0.01 * (random.randint(-90, 90))
        delta = math.ceil(min(exon_len, intron_len) * splice_perc)

        new_transcript_id = self._gv.duplicate_transcript(transcript_id)
        lh.debug(f"{new_transcript_id}: A3P {exon_id}, {delta}")
        new_transcript = self._gv.get_transcript(new_transcript_id)
        new_transcript.get_nth_exon(exon_id).start += delta
        self._gv.transcript_sort_exons(new_transcript_id)
        return new_transcript_id

    def perform_alternative_5p_splicing(self, gene: Gene):
        self.general_perform_wrapper(gene, self.core_perform_alternative_5p_splicing)

    def perform_alternative_splicing(self, gene: Gene) -> None:
        func: Callable[[Gene], None] = random.choices(
            population=(
                self.perform_alternative_3p_splicing,
                self.perform_intron_retention,
                self.perform_alternative_5p_splicing,
                self.perform_exon_skipping,
                lambda _gene: map(
                    lambda f: f(_gene), [
                        lambda x: self.perform_alternative_5p_splicing(x),
                        lambda x: self.perform_alternative_3p_splicing(x)
                    ]
                ),
                lambda _gene: map(
                    lambda f: f(_gene), [
                        lambda x: self.perform_exon_skipping(x),
                        lambda x: self.perform_alternative_3p_splicing(x)
                    ]
                ),
                lambda _gene: map(
                    lambda f: f(_gene), [
                        lambda x: self.perform_intron_retention(x),
                        lambda x: self.perform_alternative_3p_splicing(x)
                    ]
                )
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
        func(gene)

    def try_generate_n_isoform_for_a_gene(self, gene: Gene, n: int):
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
            while gene.number_of_transcripts < n:
                try:
                    self.perform_alternative_splicing(gene)
                except ImpossibleToGenerateASEventError:
                    number_of_fail += 1
                if number_of_fail > 2 * n:
                    raise ImpossibleToGenerateASEventError("Generation FAILED!")
                elif gene.number_of_transcripts == n:
                    return

    def run(self, mu: Union[int, float]):
        gene_ids_to_del: List[str] = []
        for gene in tqdm(
                iterable=list(self._gv.iter_genes()),
                desc="Generating isoforms...",
                total=self._gv.number_of_genes
        ):
            n = max(0, min(int(random.gauss(mu, mu / 2)), mu * 2))
            try:
                self.try_generate_n_isoform_for_a_gene(gene, n)
            except ImpossibleToGenerateASEventError:
                gene_ids_to_del.append(gene.gene_id)

        lh.info(f"Will remove {len(gene_ids_to_del)} genes out of {self._gv.number_of_genes}")
        for gene_id in gene_ids_to_del:
            self._gv.del_gene(gene_id)
        lh.info(f"Will remove genes FIN")
        self._gv.standardize()

    def to_file(self, filename: str):
        with open(filename + ".n_of_transcript_in_a_gene.tsv", "w") as tn:
            tn.write('gene' + '\t' + 'number_of_transcripts' + "\n")
            for gene in self._gv.iter_genes():
                tn.write(str(gene) + '\t' + str(gene.number_of_transcripts) + "\n")
        self._gv.to_file(filename)


if __name__ == '__main__':
    for i in range(1, 10, 2):
        gv = GeneViewFactory.from_file("/media/yuzj/BUP/iter_4/ce11.ncbiRefSeq.gtf")
        asm = ASManipulator(gv)
        asm.run(i)
        asm.to_file(f"{i}.gtf")
