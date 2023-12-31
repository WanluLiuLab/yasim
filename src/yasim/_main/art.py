"""
art.py -- LLRG adapter for ART, a NGS DNA-Seq simulator

.. versionadded:: 3.1.5
"""

__all__ = ("main", "create_parser")

import argparse

from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.typing_importer import List
from yasim.helper import llrg
from yasim.helper.rna_seq import bulk_rna_seq_frontend
from yasim.llrg_adapter import art


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim art",
        description=__doc__.splitlines()[1],
    )
    parser = llrg.patch_frontend_parser_public(
        parser,
        llrg_name="art",
        default_llrg_executable_name="art_illumina",
    )
    parser = llrg.patch_frontend_parser_bulk_rna_seq(parser)
    parser = art.patch_frontend_parser(parser)
    return parser


def main(args: List[str]) -> int:
    args, other_args = create_parser().parse_known_args(args)

    return bulk_rna_seq_frontend(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="art_illumina" if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "other_args": other_args,
            "sequencer_name": args.sequencer_name,
            "read_length": args.read_length,
            "pair_end_fragment_length_mean": args.pair_end_fragment_length_mean,
            "pair_end_fragment_length_std": args.pair_end_fragment_length_std,
            "is_pair_end": args.is_pair_end,
            "preserve_intermediate_files": args.preserve_intermediate_files,
        },
        assembler_args={},
        adapter_class=art.ArtAdapter,
        is_pair_end=args.is_pair_end,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=args.not_perform_assemble,
    )
