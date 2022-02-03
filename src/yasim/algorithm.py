import os.path
import random
from typing import Dict

import commonutils.parallel_helper
from bioutils.datastructure.gene import GeneView
from bioutils.io.fasta import FastaView
from commonutils import ioctl
from commonutils.logger import get_logger
from commonutils.tqdm_importer import tqdm

logger = get_logger(__name__)


def sample_exon(
        input_gtf_filename: str,
        output_gtf_filename: str,
        fasta_handler: FastaView
) -> GeneView:
    logger.info("Creating Geneview")
    gv = GeneView(input_gtf_filename)
    logger.info("Creating Geneview FIN")
    logger.info(f"Loaded {len(gv.genes)} genes with {len(gv.transcripts)} transcript")
    transcript_to_del = []
    for k, v in tqdm(iterable=gv.transcripts.items(), desc="Sampling Exons..."):
        indices = random.sample(range(len(v.exons)), int(len(v.exons)*0.75))
        v.exons = [v.exons[i] for i in sorted(indices)]
        if len(v.transcribe(fasta_handler)) >= 250:
            pass
        else:
            transcript_to_del.append(k)
    for transcript in tqdm(iterable=transcript_to_del, desc="Deleting unwanted transcripts"):
        gv.del_transcript(transcript)
    logger.info(f"Remaining {len(gv.genes)} genes with {len(gv.transcripts)} transcript")
    gv.to_file(output_gtf_filename)
    return gv


def transcribe(
        gv: GeneView,
        output_fasta: str,
        fasta_handler: FastaView
):
    output_dirname = output_fasta+".d"
    ioctl.mkdir_p(output_dirname)
    with ioctl.get_writer(output_fasta) as writer:
        for k, v in tqdm(iterable=gv.transcripts.items(), desc="Transcribing GTF..."):
            fa_name = k
            fa_value = v.transcribe(fasta_handler)
            fa_str = f">{fa_name}\n{fa_value}\n"
            writer.write(fa_str)
            with ioctl.get_writer(os.path.join(output_dirname, fa_name+".fa")) as split_writer:
                split_writer.write(fa_str)

def simulate_dge(
        gv: GeneView,
        output_tsv:str,
        max_expression:int=100
) -> Dict[str, int]:
    gene_names = gv.genes.keys()
    depth={}
    with ioctl.get_writer(output_tsv) as writer:
        writer.write(f"gene_name\tdepth\n")
        for gene_name in gene_names:
            d = int(random.uniform(0, max_expression))
            writer.write(f"{gene_name}\t{d}\n")
            depth[gene_name] = d
    return depth

def simulate(
        depth: Dict[str, int],
        output_fasta_dir: str,

):
    simulating_pool = commonutils.parallel_helper.ParallelJobQueue(pool_name="Simulating jobs")
    for gene_name in depth.keys():
        if not ioctl.file_exists(os.path.join(output_fasta_dir, gene_name+".fa")):
            pass
        simulator
    simulating_pool.start()
    simulating_pool.join()

if __name__ == "__main__":
    logger.info("Creating FastaView")
    fv = FastaView("hg38.fa")
    gv = sample_exon("hg38.knownGene.gtf", "hg38.knownGene.shuffled.gtf", fv)
    transcribe(gv, "hg38.transcripts.fa", fv)
