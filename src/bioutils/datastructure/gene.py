import os
from typing import Dict

from bioutils.datastructure import Gene, Transcript, SimpleExon
from bioutils.io import gtf
from bioutils.io.gtf import GtfRecord
from commonutils import ioctl
from commonutils import pickle_with_tqdm
from commonutils.logger import get_logger
from commonutils.tqdm_importer import tqdm

lh = get_logger(__name__)


class GeneView:
    gene_filename: str
    genes: Dict[str, Gene]
    transcripts: Dict[str, Transcript]

    def __init__(
            self,
            gene_filename: str,
            file_type: str = "gtf"
    ):
        index_filename = gene_filename + ".gvpkl.xz"
        if ioctl.file_exists(index_filename):
            if os.path.getmtime(index_filename) - os.path.getmtime(gene_filename) > 0:
                try:
                    (self.genes, self.transcripts) = pickle_with_tqdm.load(index_filename)
                    return
                except Exception:
                    lh.error("Gene index broken, will rebuild one.")
            else:
                lh.log("Gene index too old, will rebuild one.")
        if file_type == "gtf":
            if ioctl.ensure_input_existence(gene_filename):
                self.gene_filename = gene_filename
                self.read_gtf()
        else:
            pass
        pickle_with_tqdm.dump((self.genes, self.transcripts), index_filename)

    def read_gtf(self):
        def add_gene(_gtf_record: GtfRecord) -> str:
            gene_id = _gtf_record.attribute['gene_id']
            if gene_id not in self.genes.keys():
                self.genes[gene_id] = Gene.from_gtf_record(_gtf_record)
            return gene_id

        def add_transcript(_gtf_record: GtfRecord) -> str:
            gene_id = add_gene(_gtf_record)

            # Try add transcript
            transcript_id = _gtf_record.attribute['transcript_id']
            if not transcript_id in self.transcripts.keys():
                self.transcripts[transcript_id] = Transcript.from_gtf_record(_gtf_record)
            transcript = self.transcripts[transcript_id]
            self.genes[gene_id].transcripts[transcript_id] = transcript
            return transcript_id

        def add_exon(_gtf_record: GtfRecord):
            transcript_id = add_transcript(_gtf_record)

            # Try add exon
            exon = SimpleExon.from_gtf_record(_gtf_record)

            # Add exon to transcript
            self.transcripts[transcript_id].exons.append(exon)

        self.genes = {}
        self.transcripts = {}
        for gtf_record in gtf.GtfIterator(self.gene_filename):
            if "gene_id" not in gtf_record.attribute.keys() or "transcript_id" not in gtf_record.attribute.keys():
                continue
            if gtf_record.feature == "gene":
                add_gene(gtf_record)
            elif gtf_record.feature == "transcript":
                add_transcript(gtf_record)
            elif gtf_record.feature == "exon":
                add_exon(gtf_record)

    def del_gene(self, name: str):
        if name in self.genes.keys():
            for transcript in self.genes[name].transcripts.keys():
                self.del_transcript(transcript)
            self.genes.pop(name)

    def del_transcript(self, name: str):
        if name in self.transcripts.keys():
            gene_name = self.transcripts[name].gene_name
            if gene_name is not None:
                self.genes[gene_name].transcripts.pop(name)
            if len(self.genes[gene_name].transcripts) == 0:
                self.genes.pop(gene_name)
            self.transcripts.pop(name)

    def to_file(self, output_filename: str):
        with ioctl.get_writer(output_filename) as writer:
            for gk, gv in tqdm(iterable=self.genes.items(), desc=f"Writing to {output_filename}"):
                if gv.data is not None and gv.data.feature == "gene":
                    writer.write(str(gv.data) + "\n")
                for tk, tv in gv.transcripts.items():
                    if tv.data is not None and tv.data.feature == "transcript":
                        writer.write(str(tv.data) + "\n")
                    for exon in tv.exons:
                        writer.write(str(exon.data) + "\n")
