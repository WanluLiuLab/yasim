"""
assemble.py -- Assemble unassembled outputs.

.. versionadded:: 3.1.5
"""

__all__ = ("main", "create_parser")

import argparse
import os
import time

from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Type
from yasim.helper import llrg, depth_io
from yasim.helper.llrg import AssemblePairEnd, AssembleSingleEnd, AssemblerType

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim assemble",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-F",
        "--fastas",
        required=True,
        help="Directory of transcribed cDNA sequences in FASTA format from `transcribe` step",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--simulator_name",
        required=True,
        help="Custom simulator name. Used in FASTQ tags.",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-d",
        "--depth",
        required=True,
        help="Path to input Isoform-Level Depth TSV generated by `generate_depth_v2` or `generate_isoform_depth` step",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output transcript prefix. Should be prefix path to {out}.d directory.",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-i",
        "--input_fastq_dir",
        required=True,
        help="Input transcript prefix. Should be prefix path to {out}.d directory.",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--is_pair_end",
        required=False,
        help="Whether to use Pair End (PE) Simulation",
        action="store_true",
    )
    parser = llrg.patch_frontend_parser_tgs(parser)
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    assembler_class: Type[AssemblerType]
    if args.is_pair_end:
        assembler_class = AssemblePairEnd
    else:
        assembler_class = AssembleSingleEnd
    depth_data = depth_io.read_depth(args.depth)
    assembler = assembler_class(
        depth_data=depth_data,
        output_fastq_prefix=os.path.abspath(args.out),
        input_fastq_dir=os.path.abspath(args.input_fastq_dir),
        simulator_name=args.simulator_name,
        input_transcriptome_fasta_dir=os.path.abspath(args.fastas),
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
    )
    assembler.start()
    for transcript_id in depth_data.keys():
        assembler.add_transcript_id(transcript_id)

    assembler.terminate()
    while assembler.is_alive():
        _lh.info("ASSEMB: %s -- PENDING: %d", args.out, assembler.n_pending)
        time.sleep(1.0)
    assembler.join()
