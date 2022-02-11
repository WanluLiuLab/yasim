import argparse
import os.path
from typing import List

from bioutils.datastructure import GeneView
from bioutils.io.fasta import FastaView
from commonutils import ioctl
from commonutils.logger import get_logger
from commonutils.tqdm_importer import tqdm

logger = get_logger(__name__)

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="Reference genome, in FASTA format", nargs='?',
                        type=str, action='store')
    parser.add_argument('-g', '--gtf', required=True, help="Input GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Name of Output FASTA", nargs='?',
                        type=str, action='store')
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)


def transcribe(
        gv: GeneView,
        output_fasta: str,
        fasta_handler: FastaView
):
    ioctl.ensure_output_existence(output_fasta)
    with ioctl.get_writer(os.path.join(output_fasta)) as writer:
        for k, v in tqdm(iterable=gv.transcripts.items(), desc="Transcribing GTF..."):
            fa_name = k
            fa_value = v.cdna_sequence(fasta_handler)
            fa_str = f">{fa_name}\n{fa_value}\n"
            writer.write(fa_str)


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneView._from_gtf(args.gtf)
    fv = FastaView(args.fasta)
    transcribe(gv, args.out, fv)
