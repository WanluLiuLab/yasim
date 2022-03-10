from __future__ import annotations

import math
import os
import time
from abc import abstractmethod
from typing import Optional, Dict, Iterator, Union, Type

from bioutils.datastructure.gene_view_proxy import Gene, Transcript, Exon
from bioutils.io import get_file_type_from_suffix
from bioutils.io.feature import GtfIterator, GtfWriter, Gff3Iterator
from bioutils.typing.feature import GtfRecord, Feature, Gff3Record
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.file_system import file_exists
from commonutils.stdlib_helper import pickle_helper
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

GVPKL_VERSION = 0.3
"""Current version of GVPKL standard."""


class BaseGeneView:
    """
    Base class of all GeneViews.
    You may inherit from this class if you wish to add your own file format.
    """

    genes: Dict[str, Gene]
    """
    Stores a mapping of gene_id -> Gene proxy object.
    """

    transcripts: Dict[str, Transcript]
    """
    Stores a mapping of transcript_id -> Transcript proxy object.
    """

    def __init__(self):
        self.genes = {}
        self.transcripts = {}

    @classmethod
    @abstractmethod
    def _from_own_filetype(cls, filename: str):
        """
        Generate index de novo.
        """
        pass

    @classmethod
    def from_file(cls,
                  filename: str,
                  not_save_index: bool = False,
                  **kwargs):
        """
        Load GeneView from file. This function does follow things:

        * Try to find an index with suffix ".gvpkl.gz" which is newer than the file.
        * If failed, build and index from file by calling :py:func:`_from_own_filetype` and save it.
            set ``not_save_index`` to ``True`` to prevent this behaviour.
        """
        index_filename = filename + ".gvpkl.xz"
        if file_exists(index_filename) and \
                os.path.getmtime(index_filename) - os.path.getmtime(filename) > 0:
            try:
                return cls._from_gvpkl(index_filename)
            except Exception:
                lh.error("Gene index broken or too old, will rebuild one.")
        new_instance = cls._from_own_filetype(filename)
        if not not_save_index:
            new_instance.to_gvpkl(index_filename)
        return new_instance

    @classmethod
    @abstractmethod
    def from_iterator(cls, iterator: Iterator[Feature]):
        """
        Build GeneView from an iterator of Feature.
        """
        pass

    @classmethod
    def _from_gvpkl(cls, index_filename: str):
        """
        Read GVPKL index files
        """
        new_instance = cls()
        (version, new_instance.genes, new_instance.transcripts) = pickle_helper.load(index_filename)
        if version != GVPKL_VERSION:
            raise ValueError("Version mismatch")
        return new_instance

    def to_gvpkl(self, index_filename: str):
        """
        Save current GeneView object to GVPKL.
        """
        lh.info("Pickling to gvpkl...")
        pickle_helper.dump((GVPKL_VERSION, self.genes, self.transcripts), index_filename)

    def standardize(self):
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
        self._standardize_transcripts()
        self._standardize_genes()

    def _standardize_transcripts(self):
        transcript_id_to_del = []
        """Transcript to be deleted for reason like no exons."""
        for transcript in tqdm(iterable=self.transcripts.values(), desc="Standardizing transcripts"):
            if len(transcript.exons) == 0:
                transcript_id_to_del.append(transcript.transcript_id)
                continue
            transcript.sort_exons()

            if transcript.feature != "transcript":
                transcript.copy_data()
                exon_s_min = math.inf
                exon_e_max = - math.inf
                for exon in transcript.exons:
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
        for gene in tqdm(iterable=self.genes.values(), desc="Standardizing genes"):
            if len(gene.transcripts) == 0:
                gene_id_to_del.append(gene.gene_id)
                continue
            if gene.feature != "gene":
                gene.copy_data()
                transcript_s_min = math.inf
                transcript_e_max = - math.inf
                for transcripts in gene.transcripts.values():
                    transcript_s_min = min(transcript_s_min, transcripts.start)
                    transcript_e_max = max(transcript_e_max, transcripts.end)
                gene.start = transcript_s_min
                gene.end = transcript_e_max
                gene.feature = "gene"
        for gene_id in gene_id_to_del:
            self.del_gene(gene_id)

    def get_iterator(self) -> Iterator[Feature]:
        """
        Get iterator for Gene-Transcript-Exon Three-Tier Structure.
        """
        for gene in self.genes.values():
            if gene.feature == "gene":
                yield gene.get_data()
            for transcript in gene.transcripts.values():
                if transcript.feature == "transcript":
                    yield transcript.get_data()
                for exon in transcript.exons:
                    yield exon.get_data()

    def del_gene(self, gene_id: str):
        """
        Remove a gene.
        """
        if gene_id in self.genes.keys():
            for transcript_id in self.genes[gene_id].transcripts.keys():
                self.del_transcript(transcript_id)
            self.genes.pop(gene_id)

    def del_transcript(self, transcript_id: str):
        """
        Remove a transcript.
        If this is the last transcript of a gene, the gene will be removed as well.
        """
        if transcript_id in self.transcripts.keys():
            gene_id = self.transcripts[transcript_id].gene_id
            if gene_id is not None:
                self.genes[gene_id].transcripts.pop(transcript_id)
            if len(self.genes[gene_id].transcripts) == 0:
                self.genes.pop(gene_id)
            self.transcripts.pop(transcript_id)

    @abstractmethod
    def to_file(self, output_filename: str):
        """
        Write GeneView to corresponding file.
        """
        pass

    def __len__(self) -> int:
        """
        Get number of records inside.
        """
        return len(list(self.get_iterator()))


class _GtfGeneView(BaseGeneView):

    @classmethod
    @abstractmethod
    def from_iterator(cls, iterator: Union[Iterator[GtfRecord], GtfIterator]):
        def register_gene(_new_instance: _GtfGeneView, record: GtfRecord):
            gene_id = record.attribute['gene_id']
            if not record.attribute['gene_id'] in _new_instance.genes.keys():
                _new_instance.genes[gene_id] = Gene.from_feature(record)

        def register_transcript(_new_instance: _GtfGeneView, record: GtfRecord):
            gene_id = record.attribute['gene_id']
            transcript_id = record.attribute['transcript_id']
            transcript = Transcript.from_feature(record)
            if gene_id not in _new_instance.genes.keys():
                register_gene(_new_instance, record)
            if transcript_id not in _new_instance.genes[gene_id].transcripts.keys():
                _new_instance.genes[gene_id].transcripts[transcript_id] = transcript
            if transcript_id not in _new_instance.transcripts.keys():
                _new_instance.transcripts[transcript_id] = transcript

        def register_exon(_new_instance: _GtfGeneView, record: GtfRecord):
            transcript_id = record.attribute['transcript_id']
            if not record.attribute['transcript_id'] in _new_instance.transcripts.keys():
                register_transcript(_new_instance, record)
            _new_instance.transcripts[transcript_id].exons.append(Exon.from_feature(record))

        new_instance = cls()

        for gtf_record in iterator:
            if gtf_record.feature == "gene":
                register_gene(new_instance, gtf_record)
            elif gtf_record.feature == 'transcript':
                register_transcript(new_instance, gtf_record)
            elif gtf_record.feature == 'exon':
                register_exon(new_instance, gtf_record)
        return new_instance

    @classmethod
    def _from_own_filetype(cls, filename: str):
        return cls.from_iterator(GtfIterator(filename))

    def to_file(self, output_filename: str):
        GtfWriter.write_iterator(
            self.get_iterator(),
            output_filename,
            [f'created by Geneview at {time.asctime()}']
        )


class _Gff3GeneView(BaseGeneView):
    pass
    # raise NotImplementedError


class GeneView(BaseGeneView):
    """
    GeneView Factory that creates GeneView according to different file types.
    """
    @classmethod
    def from_file(cls, filename: str, file_type:Optional[str]=None, **kwargs) -> GeneView:
        if file_type is None:
            file_type = get_file_type_from_suffix(filename)
        if file_type == "GTF":
            return _GtfGeneView.from_file(filename, **kwargs)
        elif file_type == "GFF3":
            return _Gff3GeneView.from_file(filename, **kwargs)
        else:
            raise ValueError(f"Unknown file type {file_type} for {filename}")

    @classmethod
    def from_iterator(cls, iterator: Iterator[Feature], record_type: Optional[Type] = None):
        if record_type == GtfRecord or isinstance(iterator, GtfIterator):
            return _GtfGeneView.from_iterator(iterator)
        elif record_type == Gff3Record or isinstance(iterator, Gff3Iterator):
            return _Gff3GeneView.from_iterator(iterator)
        else:
            raise ValueError(f"Unknown iterator type {type(iterator)}")
