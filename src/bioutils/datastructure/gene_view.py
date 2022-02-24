__version__ = 0.2

import os
import time
from typing import Optional, Dict, Iterator

import commonutils.io.file_system
from bioutils.datastructure.gene_typing import Gene, Transcript, Exon
from bioutils.datastructure.gff_gtf_record import GtfRecord, Gff3Record
from bioutils.io import get_file_type_from_suffix
from bioutils.io.feature import Gff3Tree, GtfIterator, GtfWriter, Gff3Writer
from commonutils import pickle_with_tqdm
from commonutils.logger import get_logger

lh = get_logger(__name__)


class GeneView:
    genes: Dict[str, Gene]
    transcripts: Dict[str, Transcript]

    @classmethod
    def from_file(cls, filename: str, file_type: Optional[str] = None):
        index_filename = filename + ".gvpkl.xz"
        if commonutils.io.file_system.file_exists(index_filename) and \
                os.path.getmtime(index_filename) - os.path.getmtime(filename) > 0:
            try:
                return cls._from_gvpkl(index_filename)
            except Exception:
                lh.error("Gene index broken or too old, will rebuild one.")
        if file_type is None:
            file_type = get_file_type_from_suffix(filename)
        if file_type == "GTF":
            return cls._from_gtf(filename)
        elif file_type == "GFF3":
            return cls._from_gff3(filename)
        else:
            raise ValueError(f"Unknown file type {file_type} for {filename}")

    def __init__(self):
        self.genes = {}
        self.transcripts = {}

    @classmethod
    def _from_gvpkl(cls, index_filename: str):
        new_instance = cls()
        (version, new_instance.genes, new_instance.transcripts) = pickle_with_tqdm.load(index_filename)
        if version != __version__:
            raise ValueError("Version mismatch")
        return new_instance

    def to_gvpkl(self, index_filename: str):
        lh.info("Pickling to gvpkl...")
        pickle_with_tqdm.dump((__version__, self.genes, self.transcripts), index_filename)

    @classmethod
    def _from_gtf(cls, filename: str):
        def register_gene(_new_instance: GeneView, record: GtfRecord):
            gene_id = record.attribute['gene_id']
            if not record.attribute['gene_id'] in _new_instance.genes.keys():
                _new_instance.genes[gene_id] = Gene.from_gtf_record(record)

        def register_transcript(_new_instance: GeneView, record: GtfRecord):
            gene_id = record.attribute['gene_id']
            transcript_id = record.attribute['transcript_id']
            transcript = Transcript.from_gtf_record(record)
            if gene_id not in _new_instance.genes.keys():
                register_gene(_new_instance, record)
            if transcript_id not in _new_instance.genes[gene_id].transcripts.keys():
                _new_instance.genes[gene_id].transcripts[transcript_id] = transcript
            if transcript_id not in _new_instance.transcripts.keys():
                _new_instance.transcripts[transcript_id] = transcript

        def register_exon(_new_instance: GeneView, record: GtfRecord):
            transcript_id = record.attribute['transcript_id']
            if not record.attribute['transcript_id'] in _new_instance.transcripts.keys():
                register_transcript(_new_instance, record)
            _new_instance.transcripts[transcript_id].exons.append(Exon.from_gtf_record(record))

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

    def standardize(self):
        return
        self.standardize_transcripts()
        self.standardize_genes()

    def standardize_transcripts(self):
        raise NotImplementedError  # TODO: Implement

    def standardize_genes(self):
        raise NotImplementedError  # TODO: Implement

    def get_transcript_iterator(self) -> Iterator[Transcript]:
        for transcript in self.transcripts.values():
            yield transcript

    @classmethod
    def _from_gff2(cls, filename: str):
        return cls._from_gtf(filename)

    @classmethod
    def _from_gff3(cls, filename: str):
        raise NotImplementedError  # See https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

        def gff3_bfs(_new_instance: GeneView, root_id: str):
            # TODO: new_instance()

            for child_id in gff3_tree.get_child_ids(root_id):
                gff3_bfs(new_instance, child_id)

        new_instance = cls()

        gff3_tree = Gff3Tree(filename)
        for root_id in gff3_tree.get_toplevel_ids():
            gff3_bfs(new_instance, root_id)

    def get_gtf_iterator(self) -> Iterator[GtfRecord]:
        for gene in self.genes.values():
            yield gene.to_gtf_record()
            for transcript in gene.transcripts.values():
                yield transcript.to_gtf_record()
                for exon in transcript.exons:
                    yield exon.to_gtf_record()

    def to_gtf(self, output_filename: str):
        GtfWriter.write_iterator(
            self.get_gtf_iterator(),
            output_filename,
            [f'created by Geneview at {time.asctime()}']
        )

    def get_gff3_iterator(self) -> Iterator[Gff3Record]:
        raise NotImplementedError  # TODO: Implement

    def to_gff3(self, output_filename: str):
        Gff3Writer.write_iterator(
            self.get_gff3_iterator(),
            output_filename,
            [f'created by Geneview at {time.asctime()}']
        )

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
