import argparse
from typing import List

from yasim.helper.rna_seq import bulk_rna_seq_frontend
from yasim.helper import llrg
from yasim.llrg_adapter import pbsim2


def main(args: List[str]):
    parser = argparse.ArgumentParser()
    parser = llrg.patch_frontend_parser_public(
        parser,
        llrg_name="pbsim2",
        default_llrg_executable_name="pbsim2"
    )
    parser = llrg.patch_frontend_parser_bulk_rna_seq(parser)
    parser = llrg.patch_frontend_parser_tgs(parser)
    parser = pbsim2.patch_frontend_parser(parser)
    args, other_args = parser.parse_known_args(args)
    return bulk_rna_seq_frontend(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        depth_file_path=args.depth,
        jobs=args.jobs,
        simulator_name="_".join(
            ("pbsim2", args.hmm_model, "clr")
        ) if args.simulator_name is None else args.simulator_name,
        adapter_args={
            "hmm_model": args.hmm_model,
            "other_args": other_args
        },
        assembler_args={
            "truncate_ratio_3p": args.truncate_ratio_3p,
            "truncate_ratio_5p": args.truncate_ratio_5p
        },
        adapter_class=pbsim2.Pbsim2Adapter,
        is_pair_end=False,
        llrg_executable_path=args.llrg_executable_path,
        not_perform_assemble=False
    )
