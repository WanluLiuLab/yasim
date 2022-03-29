import argparse
import os
from typing import List

from bioutils.algorithm.sequence import get_gc_percent
from bioutils.datastructure.fasta_view import FastaView
from bioutils.datastructure.gene_view import GeneViewFactory, GeneViewType
from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer
from commonutils.stdlib_helper.logger_helper import get_logger

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Name of Output FASTA", nargs='?',
                        type=str, action='store')
    return parser.parse_args(args)


def transcribe(
        gv: GeneViewType,
        output_fasta: str,
        fv: FastaView
):
    intermediate_fasta_dir = output_fasta + ".d"
    shell_utils.mkdir_p(intermediate_fasta_dir)
    with get_writer(output_fasta) as fasta_writer, \
            get_writer(output_fasta + ".stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "GENE_ID",
            "SEQNAME",
            "START",
            "END",
            "STRAND",
            "LEN",
            "GC"
        )) + "\n")
        for transcript_name, transcript_value in tqdm(iterable=gv.transcripts.items(), desc="Transcribing GTF..."):
            cdna_seq = transcript_value.cdna_sequence(sequence_func=fv.sequence)
            fa_str = f">{transcript_name}\n{cdna_seq}\n"
            fasta_writer.write(fa_str)
            stats_writer.write("\t".join((
                transcript_name,
                transcript_value.gene_id,
                transcript_value.seqname,
                str(transcript_value.start),
                str(transcript_value.end),
                transcript_value.strand,
                str(transcript_value.end - transcript_value.start),
                str(round(get_gc_percent(cdna_seq) * 100, 2))
            )) + "\n")
            transcript_output_fasta = os.path.join(intermediate_fasta_dir, f"{transcript_name}.fa")
            with get_writer(transcript_output_fasta) as single_transcript_writer:
                single_transcript_writer.write(fa_str)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    fv = FastaView(args.fasta)
    transcribe(gv, args.out, fv)
