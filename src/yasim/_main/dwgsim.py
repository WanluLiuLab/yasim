import argparse
from typing import List

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.rna_seq import bulk_rna_seq_frontend
from yasim.helper import llrg
from yasim.llrg_adapter import dwgsim

logger = get_logger(__name__)


def main(args: List[str]):
    parser = argparse.ArgumentParser()
    parser = llrg.patch_frontend_parser_public(
        parser,
        llrg_name="dwgsim",
        default_llrg_executable_name="dwgsim"
    )
    parser = llrg.patch_frontend_parser_bulk_rna_seq(parser)
    args, other_args = parser.parse_known_args(args)
    return bulk_rna_seq_frontend(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="dwgsim" if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "other_args": other_args
        },
        assembler_args={},
        adapter_class=dwgsim.DwgsimAdapter,
        is_pair_end=True,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=False
    )
