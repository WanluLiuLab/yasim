import argparse
import os.path
from typing import List, Dict, Optional

from bioutils.datastructure.gene_view import GeneView
from bioutils.io.fastx import FastaView
from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer
from commonutils.stdlib_helper.logger_helper import get_logger
from yasim.main import dge

logger = get_logger(__name__)

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--depth', required=False, help="Depth generated by `dge` step", nargs='?',
                        type=str, action='store', default=None)
    parser.add_argument('-o', '--out', required=True, help="Name of Output FASTA", nargs='?',
                        type=str, action='store')
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)


def transcribe(
        gv: GeneView,
        output_fasta: str,
        fv: FastaView,
        depth: Optional[Dict[str, int]] = None
):
    with get_writer(os.path.join(output_fasta)) as writer:
        for k, v in tqdm(iterable=gv.transcripts.items(), desc="Transcribing GTF..."):
            fa_name = k
            fa_value = v.cdna_sequence(sequence_func=fv.sequence)
            fa_str = f">{fa_name}\n{fa_value}\n"
            writer.write(fa_str)
    if depth is None:
        return
    depth_cluster = dge.cluster_depth(depth)
    intermediate_fasta_dir = output_fasta + ".d"
    shell_utils.mkdir_p(intermediate_fasta_dir)
    ofv = FastaView(output_fasta)
    for transcript_depth, transcript_names in tqdm(iterable=depth_cluster.items(), desc="Transcribing DGE..."):
        transcript_output_fasta = os.path.join(intermediate_fasta_dir, f"{transcript_depth}.fa")
        try:
            ofv.subset_chr(transcript_output_fasta, transcript_names)
        except ValueError:
            pass


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneView.from_file(args.gtf)
    fv = FastaView(args.fasta)
    if args.depth is not None:
        depth = dge.read_depth(args.depth)
    else:
        depth = None
    transcribe(gv, args.out, fv, depth)
