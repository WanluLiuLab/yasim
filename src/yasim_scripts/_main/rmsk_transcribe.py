"""
TODO
"""
import argparse

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.parser.gtf import GtfIterator
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_scripts rmsk_transcribe",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-i",
        "--src_rmsk_gtf_path",
        type=str,
        required=True,
        help="File path of source RepeatMasker GTF",
    )
    parser.add_argument(
        "-f",
        "--src_fa_path",
        type=str,
        required=True,
        help="File path of reference FASTA",
    )
    parser.add_argument(
        "-o",
        "--dst_fa_path",
        type=str,
        required=True,
        help="File path of destination FASTA",
    )
    return parser


def transcribe_rmsk(src_rmsk_gtf_path: str, src_fa_path: str, dst_fa_path: str):
    with FastaWriter(dst_fa_path) as faw:
        fav = FastaViewFactory(src_fa_path)
        for record in GtfIterator(src_rmsk_gtf_path):
            faw.write(
                FastaRecord(
                    record.attribute_get("transcript_id"), fav.sequence(record.seqname, record.start0b, record.end0b)
                )
            )


def main(args: List[str]):
    argv = create_parser().parse_args(args)
    transcribe_rmsk(
        src_rmsk_gtf_path=argv.src_rmsk_gtf_path,
        src_fa_path=argv.src_fa_path,
        dst_fa_path=argv.dst_fa_path,
    )
