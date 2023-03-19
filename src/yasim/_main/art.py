__all__ = (
    "main",
)

import argparse
from typing import List

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.rna_seq import bulk_rna_seq_frontend
from yasim.helper import llrg
from yasim.llrg_adapter import art

_lh = get_logger(__name__)


def main(args: List[str]) -> int:
    parser = argparse.ArgumentParser()
    parser = llrg.patch_frontend_parser_public(
        parser,
        llrg_name="art",
        default_llrg_executable_name="art_illumina"
    )
    parser = llrg.patch_frontend_parser_bulk_rna_seq(parser)
    parser = art.patch_frontend_parser(parser)
    args, other_args = parser.parse_known_args(args)

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
            "is_pair_end": args.is_pair_end
        },
        assembler_args={},
        adapter_class=art.ArtAdapter,
        is_pair_end=args.is_pair_end,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=False
    )
