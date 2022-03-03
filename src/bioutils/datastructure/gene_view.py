__version__ = 0.3

import os
import time
from abc import abstractmethod
from typing import Optional, Dict, Iterator, Union

import commonutils.io.file_system
from bioutils.datastructure._gene_view_proxy import Gene, Transcript, Exon
from bioutils.io import get_file_type_from_suffix
from bioutils.io.feature import Gff3Tree, GtfIterator, GtfWriter, Gff3Writer
from bioutils.typing.feature import GtfRecord, Gff3Record, Feature
from commonutils.stdlib_helper import pickle_helper
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

class _BaseGeneView:
    genes: Dict[str, Gene]
    transcripts: Dict[str, Transcript]

    def __init__(self):
        self.genes = {}
        self.transcripts = {}


    def get_transcript_iterator(self) -> Iterator[Transcript]:
        for transcript in self.transcripts.values():
            yield transcript


    @classmethod
    @abstractmethod
    def _from_own_filetype(cls, filename:str):
        """
        Generate index de novo.
        """
        pass

    @classmethod
    def from_file(cls, filename:str):
        index_filename = filename + ".gvpkl.xz"
        if commonutils.io.file_system.file_exists(index_filename) and \
                os.path.getmtime(index_filename) - os.path.getmtime(filename) > 0:
            try:
                return cls._from_gvpkl(index_filename)
            except Exception:
                lh.error("Gene index broken or too old, will rebuild one.")
        return cls._from_own_filetype(filename)

    @classmethod
    def _from_gvpkl(cls, index_filename: str):
        new_instance = cls()
        (version, new_instance.genes, new_instance.transcripts) = pickle_helper.load(index_filename)
        if version != __version__:
            raise ValueError("Version mismatch")
        return new_instance

    def to_gvpkl(self, index_filename: str):
        lh.info("Pickling to gvpkl...")
        pickle_helper.dump((__version__, self.genes, self.transcripts), index_filename)


    def standardize(self):
        return
        self.standardize_transcripts()
        self.standardize_genes()

    def standardize_transcripts(self):
        raise NotImplementedError  # TODO: Implement

    def standardize_genes(self):
        raise NotImplementedError  # TODO: Implement


    def get_iterator(self) -> Iterator[Feature]:
        for gene in self.genes.values():
            yield gene._data
            for transcript in gene.transcripts.values():
                yield transcript._data
                for exon in transcript.exons:
                    yield exon._data


    def del_gene(self, gene_id: str):
        if gene_id in self.genes.keys():
            for transcript_id in self.genes[gene_id].transcripts.keys():
                self.del_transcript(transcript_id)
            self.genes.pop(gene_id)

    def del_transcript(self, transcript_id: str):
        if transcript_id in self.transcripts.keys():
            gene_id = self.transcripts[transcript_id].gene_id
            if gene_id is not None:
                self.genes[gene_id].transcripts.pop(transcript_id)
            if len(self.genes[gene_id].transcripts) == 0:
                self.genes.pop(gene_id)
            self.transcripts.pop(transcript_id)


class _GtfGeneView(_BaseGeneView):

    @classmethod
    def _from_own_filetype(cls, filename:str):
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

        for gtf_record in GtfIterator(filename):
            if gtf_record.feature == "gene":
                register_gene(new_instance, gtf_record)
            elif gtf_record.feature == 'transcript':
                register_transcript(new_instance, gtf_record)
            elif gtf_record.feature == 'exon':
                register_exon(new_instance, gtf_record)
        new_instance.standardize()
        index_filename = filename + ".gvpkl.xz"
        new_instance.to_gvpkl(index_filename)
        return new_instance


    def to_file(self, output_filename: str):
        GtfWriter.write_iterator(
            self.get_iterator(),
            output_filename,
            [f'created by Geneview at {time.asctime()}']
        )

class _Gff3GeneView(_BaseGeneView):
    pass
    # raise NotImplementedError


class GeneView(_BaseGeneView):
    @classmethod
    def from_file(cls, filename: str, file_type: Optional[str] = None) -> Union[_GtfGeneView, _Gff3GeneView]:
        if file_type is None:
            file_type = get_file_type_from_suffix(filename)
        if file_type == "GTF":
            return _GtfGeneView.from_file(filename)
        elif file_type == "GFF3":
            return _Gff3GeneView.from_file(filename)
        else:
            raise ValueError(f"Unknown file type {file_type} for {filename}")


