import argparse
from typing import List

from yasim._main import abstract_simulate
from yasim.helper.llrg import patch_frontend_parser_tgs, \
    patch_frontend_parser_public
from yasim.llrg_adapter import pbsim3 as pbsim3


def main(args: List[str]):
    parser = argparse.ArgumentParser()
    parser = patch_frontend_parser_public(
        parser,
        llrg_name="pbsim3",
        default_llrg_executable_name="pbsim3"
    )
    parser = patch_frontend_parser_tgs(parser)
    parser = pbsim3.patch_frontend_parser(parser)
    args, other_args = parser.parse_known_args(args)
    abstract_simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="_".join(
            ("pbsim3", args.hmm_model, "ccs" if args.ccs_pass > 1 else "clr")
        ) if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "hmm_model": args.hmm_model,
            "ccs_num_threads": 1 if args.ccs_pass > 1 else None,
            "samtools_executable_path": args.samtools_path if args.ccs_pass > 1 else None,
            "ccs_executable_path": args.ccs_path if args.ccs_pass > 1 else None,
            "ccs_pass": args.ccs_pass,
            "hmm_method": args.hmm_method,
            "strategy": args.strategy,
            "other_args": other_args
        }, assembler_args={
            "truncate_ratio_3p": args.truncate_ratio_3p,
            "truncate_ratio_5p": args.truncate_ratio_5p
        },
        adapter_class=pbsim3.Pbsim3Adapter,
        is_pair_end=False,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=False
    )
