import argparse
from typing import List

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim._main import abstract_simulate
from yasim.helper.llrg import patch_frontend_parser_tgs, patch_frontend_parser_public
from yasim.llrg_adapter import badread

logger = get_logger(__name__)


def main(args: List[str]):
    parser = argparse.ArgumentParser()
    parser = patch_frontend_parser_tgs(parser)
    parser = patch_frontend_parser_public(parser, "badread", "badread")
    parser = badread.patch_frontend_parser(parser)
    args, other_args = parser.parse_known_args(args)
    return abstract_simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="_".join(
            ("badread", args.model_name, "cDNA")
        ) if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "model_name": args.model_name,
            "other_args": other_args
        },
        assembler_args={
            "truncate_ratio_3p": args.truncate_ratio_3p,
            "truncate_ratio_5p": args.truncate_ratio_5p
        },
        adapter_class=badread.BadReadAdapter,
        is_pair_end=False,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=False
    )
