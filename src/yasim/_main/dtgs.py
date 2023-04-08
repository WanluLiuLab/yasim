"""
dwgsim.py -- LLRG adapter for DWGSIM, a NGS DNA-Seq simulator
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
from typing import List

from yasim.helper import llrg
from yasim.helper.rna_seq import bulk_rna_seq_frontend
from yasim.llrg_adapter import dwgsim, dtgs


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m yasim dtgs", description=__doc__.splitlines()[1])
    parser = llrg.patch_frontend_parser_public(
        parser,
        llrg_name="dtgs"
    )
    parser = llrg.patch_frontend_parser_bulk_rna_seq(parser)
    return parser


def main(args: List[str]):
    args, other_args = create_parser().parse_known_args(args)
    return bulk_rna_seq_frontend(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="dtgs" if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "other_args": other_args
        },
        assembler_args={},
        adapter_class=dtgs.DTGSAdapter,
        is_pair_end=False,
        not_perform_assemble=False
    )
