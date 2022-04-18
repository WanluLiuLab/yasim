from __future__ import annotations

import copy
import math
import os
import time
import uuid
from abc import abstractmethod, ABC
from typing import Optional, Dict, Iterator, Union, Type, Iterable

from bioutils.datastructure._gv_feature_proxy_mutator import GeneMutator, TranscriptMutator
from bioutils.datastructure.gv_feature_proxy import Gene, Transcript, Exon, BaseFeatureProxy, \
    DEFAULT_SORT_EXON_EXON_STRAND_POLICY
from bioutils.io import get_file_type_from_suffix
from bioutils.io.feature import GtfIterator, GtfWriter, Gff3Iterator
from bioutils.typing.feature import GtfRecord, FeatureType, Gff3Record
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.file_system import file_exists
from commonutils.stdlib_helper import pickle_helper
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

GVPKL_VERSION = "0.4"
"""Current version of GVPKL standard."""


class GeneViewType:
    """
    Abstract class for factory.
    """

    _genes: Dict[str, Gene]
    """
    Stores a mapping of gene_id -> Gene proxy object.
    """

    _transcripts: Dict[str, Transcript]
    """
    Stores a mapping of transcript_id -> Transcript proxy object.
    """

    @classmethod
    @abstractmethod
    def _from_own_filetype(cls, filename: str, fast: bool):
        """
        Generate index de novo.
        """
        pass

    @classmethod
    @abstractmethod
    def from_file(cls,
                  filename: str,
                  not_save_index: bool = False,
                  fast: bool = True,
                  **kwargs):
        """
        Load GeneView from file. This function does follow things:

        * Try to find an index with suffix ".gvpkl.gz" which is newer than the file.
        * If failed, build and index from file by calling :py:func:`_from_own_filetype` and save it.
            set ``not_save_index`` to ``True`` to prevent this behaviour.
        """
        pass

    @classmethod
    @abstractmethod
    def from_iterable(cls, iterable: Iterable[FeatureType], fast: bool = True):
        """
        Build GeneView from an iterable of FeatureType.
        """
        pass

    @abstractmethod
    def iter_transcripts(self):
        pass

    @abstractmethod
    def iter_transcript_ids(self):
        pass

    @abstractmethod
    def get_transcript(self, transcript_id: str) -> Transcript:
        pass

    @property
    @abstractmethod
    def number_of_transcripts(self) -> int:
        pass

    @abstractmethod
    def iter_genes(self):
        pass

    @abstractmethod
    def iter_gene_ids(self):
        pass

    @abstractmethod
    def get_gene(self, gene_id: str) -> Gene:
        pass

    @property
    @abstractmethod
    def number_of_genes(self) -> int:
        pass

    @abstractmethod
    def to_gvpkl(self, index_filename: str):
        """
        Save current GeneView object to GVPKL.
        """
        pass

    @abstractmethod
    def standardize(
            self,
            sort_exon_exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY,
            *args, **kwargs
    ):
        """
        This function standardizes GeneView into a Gene-Transcript-Exon Three-Tier Structure,
        with other functions introduced below:

        1. Transcript level
        1). Remove transcript without exons;
        2). Sort exons and re-mark exon number.
        3). If transcript is not built with feature ``transcript``,
            normalize its starting and ending position to its span length.
        2. Same things to be done for gene.
        """
        pass

    @abstractmethod
    def get_iterator(self) -> Iterator[FeatureType]:
        """
        Get iterator for Gene-Transcript-Exon Three-Tier Structure.
        """
        pass

    @abstractmethod
    def transcript_sort_exons(
            self,
            transcript_id: str,
            sort_exon_exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY
    ):
        """
        Re-index and sort exons of a transcript
        """
        pass

    @abstractmethod
    def to_file(self, output_filename: str):
        """
        Write GeneView to corresponding file.
        """
        pass

    @abstractmethod
    def del_gene(self, gene_id: str):
        """
        Remove a gene.
        """
        pass

    @abstractmethod
    def del_exon(self, transcript_id: str, exon_number: int):
        pass

    @abstractmethod
    def del_transcript(self, transcript_id: str, auto_remove_empty_gene: bool = True):
        """
        Remove a transcript.
        If this is the last transcript of a gene, the gene will be removed as well.
        """
        pass

    @abstractmethod
    def add_gene(self, gene: Gene):
        """
        Register a new gene.
        """
        pass

    @abstractmethod
    def duplicate_transcript(self, transcript_id: str) -> str:
        """
        Duplicate an existing transcript.

        :return: New transcript ID.
        """
        pass

    @abstractmethod
    def add_transcript(self, transcript: Transcript, fast: bool = False):
        """
        Register a new transcript. Will register corresponding gene if not exist.
        """
        pass

    @abstractmethod
    def add_exon(self, exon: Exon, fast: bool = False):
        """
        Register a new exon. Will register corresponding gene or transcript if not exist.
        """
        pass

    @abstractmethod
    def __len__(self) -> int:
        """
        Get number of records inside.
        """
        pass


class BaseGeneView(GeneViewType, ABC):
    """
    Base class of all GeneViews.
    You may inherit from this class if you wish to add your own file format.
    """

    def __init__(self):
        self._genes = {}
        self._transcripts = {}

    @classmethod
    def from_file(cls,
                  filename: str,
                  not_save_index: bool = False,
                  fast: bool = True,
                  **kwargs):
        index_filename = f"{filename}.{GVPKL_VERSION}.gvpkl.xz"
        if file_exists(index_filename) and \
                os.path.getmtime(index_filename) - os.path.getmtime(filename) > 0:
            try:
                return cls._from_gvpkl(index_filename)
            except Exception:
                lh.error("Gene index broken or too old, will rebuild one.")
        new_instance = cls._from_own_filetype(filename, fast=fast)
        if not not_save_index:
            new_instance.to_gvpkl(index_filename)
        return new_instance

    @classmethod
    def _from_gvpkl(cls, index_filename: str):
        new_instance = cls()
        if GVPKL_VERSION.endswith("unstable"):
            raise ValueError("Version mismatch")
        (version, new_instance._genes, new_instance._transcripts) = pickle_helper.load(index_filename)
        if version != GVPKL_VERSION:
            raise ValueError("Version mismatch")
        return new_instance

    def to_gvpkl(self, index_filename: str):
        lh.info("Pickling to gvpkl...")
        if GVPKL_VERSION.endswith("unstable"):
            return
        pickle_helper.dump((GVPKL_VERSION, self._genes, self._transcripts), index_filename)

    def iter_transcripts(self):
        return self._transcripts.values()

    def iter_transcript_ids(self):
        return self._transcripts.keys()

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self._transcripts[transcript_id]

    @property
    def number_of_transcripts(self) -> int:
        return len(self._transcripts)

    def iter_genes(self):
        return self._genes.values()

    def iter_gene_ids(self):
        return self._genes.keys()

    def get_gene(self, gene_id: str) -> Gene:
        return self._genes[gene_id]

    @property
    def number_of_genes(self) -> int:
        return len(self._genes)

    def transcript_sort_exons(
            self,
            transcript_id: str,
            sort_exon_exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY
    ):
        TranscriptMutator.sort_exons(
            self.get_transcript(transcript_id),
            exon_number_policy=sort_exon_exon_number_policy
        )

    def standardize(
            self,
            sort_exon_exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY,
            remove_transcript_without_exons: bool = True,
            rescale_inferred_transcript_from_exon_boundaries: bool = True,
            *args, **kwargs
    ):
        self._standardize_transcripts(
            sort_exon_exon_number_policy,
            remove_transcript_without_exons,
            rescale_inferred_transcript_from_exon_boundaries
        )
        self._standardize_genes()

    def _standardize_transcripts(
            self,
            sort_exon_exon_number_policy: str,
            remove_transcript_without_exons: bool,
            rescale_inferred_transcript_from_exon_boundaries: bool
    ):
        transcript_id_to_del = []
        """Transcript to be deleted for reason like no exons."""
        for transcript in tqdm(iterable=self.iter_transcripts(), desc="Standardizing transcripts"):
            if remove_transcript_without_exons and transcript.number_of_exons == 0:
                transcript_id_to_del.append(transcript.transcript_id)
                continue
            self.transcript_sort_exons(transcript.transcript_id)

            if rescale_inferred_transcript_from_exon_boundaries and \
                    transcript.feature != "transcript" and \
                    transcript.number_of_exons > 0:
                transcript.copy_data()
                exon_s_min = math.inf
                exon_e_max = - math.inf
                for exon in transcript.iter_exons():
                    exon_s_min = min(exon_s_min, exon.start)
                    exon_e_max = max(exon_e_max, exon.end)
                transcript.start = exon_s_min
                transcript.end = exon_e_max
                transcript.feature = "transcript"
        for transcript_id in transcript_id_to_del:
            self.del_transcript(transcript_id)

    def _standardize_genes(self):
        gene_id_to_del = []
        """Genes to be deleted for reason like no transcripts"""
        for gene in tqdm(iterable=self.iter_genes(), desc="Standardizing genes"):
            if gene.number_of_transcripts == 0:
                gene_id_to_del.append(gene.gene_id)
                continue
            if gene.feature != "gene":
                gene.copy_data()
                transcript_s_min = math.inf
                transcript_e_max = - math.inf
                for transcript in gene.iter_transcripts():
                    transcript_s_min = min(transcript_s_min, transcript.start)
                    transcript_e_max = max(transcript_e_max, transcript.end)
                gene.start = transcript_s_min
                gene.end = transcript_e_max
                gene.feature = "gene"
        for gene_id in gene_id_to_del:
            self.del_gene(gene_id)

    def get_iterator(self) -> Iterator[FeatureType]:
        for gene in self.iter_genes():
            if gene.feature == "gene":
                yield gene.get_data()
            for transcript in gene.iter_transcripts():
                if transcript.feature == "transcript":
                    yield transcript.get_data()
                for exon in transcript.iter_exons():
                    yield exon.get_data()

    def del_gene(self, gene_id: str):
        if gene_id in self.iter_gene_ids():
            for transcript_id in list(self.get_gene(gene_id).iter_transcript_ids()):
                self.del_transcript(transcript_id)
        else:
            raise ValueError(f"Gene ID {gene_id} not found!")

    def duplicate_transcript(self, transcript_id: str) -> str:
        transcript = self.get_transcript(transcript_id)
        new_transcript = copy.deepcopy(transcript)
        new_transcript.attribute['reference_transcript_id'] = new_transcript.transcript_id
        new_transcript.transcript_id = transcript.gene_id + str(uuid.uuid4())
        for exon in new_transcript.iter_exons():
            exon.transcript_id = new_transcript.transcript_id
        self.add_transcript(new_transcript, fast=True)
        return new_transcript.transcript_id

    def del_exon(self, transcript_id: str, exon_number: int):
        TranscriptMutator.del_exon(self.get_transcript(transcript_id), exon_number)

    def _del_transcript(self, transcript_id: str):
        self._transcripts.pop(transcript_id)

    def _del_gene(self, gene_id: str):
        self._genes.pop(gene_id)

    def del_transcript(self, transcript_id: str, auto_remove_empty_gene: bool = True):
        if transcript_id in self.iter_transcript_ids():
            gene_id = self.get_transcript(transcript_id).gene_id
            GeneMutator.del_transcript(self.get_gene(gene_id), transcript_id)
            if self.get_gene(gene_id).number_of_transcripts == 0 and auto_remove_empty_gene:
                lh.debug(f"Automatically remove empty gene {gene_id}")
                self._del_gene(gene_id)
            self._del_transcript(transcript_id)
        else:
            raise ValueError(f"Transcript ID {transcript_id} not found!")

    def add_gene(self, gene: Gene):
        gene_id = gene.gene_id
        if gene_id not in self.iter_gene_ids() or self.get_gene(gene_id).feature != "gene":
            self._add_gene(gene)
        if gene.feature != "gene":
            lh.warn(f"Gene {gene_id} is inferred from feature {gene.feature}")

    def _add_gene(self, gene: Gene):
        self._genes[gene.gene_id] = gene

    def _add_transcript(self, transcript: Transcript):
        self._transcripts[transcript.transcript_id] = transcript

    def add_transcript(
            self,
            transcript: Transcript,
            fast: bool = False
    ):
        gene_id = transcript.gene_id
        transcript_id = transcript.transcript_id
        if gene_id not in self.iter_gene_ids():
            self.add_gene(BaseFeatureProxy.duplicate_cast(transcript, Gene))
        gene = self.get_gene(gene_id)
        if transcript_id not in gene.iter_transcript_ids() or gene.get_transcript(
                transcript_id).feature != "transcript":
            if fast:
                GeneMutator.fast_add_transcript(gene, transcript)
            else:
                GeneMutator.add_transcript(gene, transcript)
        if transcript_id not in self.iter_transcript_ids() or \
                self.get_transcript(transcript_id).feature != "transcript":
            self._add_transcript(transcript)
        if transcript.feature != "transcript":
            lh.warn(f"Transcript {transcript_id} is inferred from feature {transcript.feature}")

    def add_exon(self, exon: Exon, fast: bool = False):
        transcript_id = exon.transcript_id
        if transcript_id not in self.iter_transcript_ids():
            self.add_transcript(BaseFeatureProxy.duplicate_cast(exon, Transcript), fast=fast)
        if fast:
            TranscriptMutator.fast_add_exon(self.get_transcript(transcript_id), exon)
        else:
            TranscriptMutator.add_exon(self.get_transcript(transcript_id), exon)

    def __len__(self) -> int:
        return len(list(self.get_iterator()))


class _GtfGeneView(BaseGeneView):

    @classmethod
    def from_iterable(cls, iterable: Union[Iterable[GtfRecord], GtfIterator], fast: bool = True):
        new_instance = cls()
        for gtf_record in iterable:
            if gtf_record.feature == "gene":
                new_instance.add_gene(Gene.from_feature(gtf_record))
            elif gtf_record.feature == 'transcript':
                new_instance.add_transcript(Transcript.from_feature(gtf_record), fast=fast)
            elif gtf_record.feature == 'exon':
                new_instance.add_exon(Exon.from_feature(gtf_record), fast=fast)
        return new_instance

    @classmethod
    def _from_own_filetype(cls, filename: str, fast: bool):
        return cls.from_iterable(GtfIterator(filename), fast=fast)

    def to_file(self, output_filename: str):
        GtfWriter.write_iterator(
            self.get_iterator(),
            output_filename,
            [f'created by Geneview at {time.asctime()}']
        )


class _Gff3GeneView(BaseGeneView):
    pass


class GeneViewFactory:
    """
    GeneView Factory that creates GeneView according to different file types.
    """

    @classmethod
    def from_file(cls, filename: str, file_type: Optional[str] = None, **kwargs) -> GeneViewType:
        if file_type is None:
            file_type = get_file_type_from_suffix(filename)
        if file_type == "GTF":
            return _GtfGeneView.from_file(filename, **kwargs)
        elif file_type == "GFF3":
            return _Gff3GeneView.from_file(filename, **kwargs)
        else:
            raise ValueError(f"Unknown file type {file_type} for {filename}")

    @classmethod
    def from_iterable(
            cls,
            iterable: Iterable[FeatureType],
            record_type: Optional[Type] = None,
            fast: bool = True
    ) -> GeneViewType:
        if record_type is None:
            record_type = type(iterable)
        if record_type == GtfRecord or isinstance(iterable, GtfIterator):
            return _GtfGeneView.from_iterable(iterable, fast=fast)
        elif record_type == Gff3Record or isinstance(iterable, Gff3Iterator):
            return _Gff3GeneView.from_iterable(iterable, fast=fast)
        else:
            raise ValueError(f"Unknown iterable type {type(iterable)}")
