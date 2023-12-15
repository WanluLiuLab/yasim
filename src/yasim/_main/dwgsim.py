"""
dwgsim.py -- LLRG adapter for DWGSIM, a NGS DNA-Seq simulator

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
from yasim.llrg_adapter import dwgsim


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim dwgsim",
        description=__doc__.splitlines()[1],
    )
    parser = llrg.patch_frontend_parser_public(
        parser,
        llrg_name="dwgsim",
        default_llrg_executable_name="dwgsim",
    )
    parser = llrg.patch_frontend_parser_bulk_rna_seq(parser)
    parser = dwgsim.patch_frontend_parser(parser)
    return parser


def main(args: List[str]):
    args, other_args = create_parser().parse_known_args(args)
    return bulk_rna_seq_frontend(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="dwgsim" if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "other_args": other_args,
            "preserve_intermediate_files": args.preserve_intermediate_files,
        },
        assembler_args={},
        adapter_class=dwgsim.DwgsimAdapter,
        is_pair_end=True,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=args.not_perform_assemble,
    )
