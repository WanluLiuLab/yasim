import argparse
from typing import List, Tuple, Optional

from yasim._main import abstract_simulate
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_tgs, \
    enhanced_which, patch_frontend_parser_public
from yasim.llrg_adapter import pbsim3 as pbsim3


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser = patch_frontend_parser_public(
        parser,
        llrg_name="pbsim3",
        default_llrg_executable_name="pbsim3"
    )
    parser = patch_frontend_parser_tgs(parser)
    parser = pbsim3.patch_frontend_parser(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        hmm_model: str,
        llrg_executable_path: str,
        depth: DepthType,
        jobs: int,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float,
        hmm_method: str,
        samtools_path: str,
        ccs_path: str,
        ccs_pass: int,
        strategy: str,
        simulator_name: Optional[str],
        other_args: List[str]
):
    abstract_simulate(
        transcriptome_fasta_dir=transcriptome_fasta_dir,
        output_fastq_prefix=output_fastq_prefix,
        depth=depth,
        jobs=jobs,
        simulator_name="_".join(
            ("pbsim3", hmm_model, "ccs" if ccs_pass > 1 else "clr")
        ) if simulator_name is None else simulator_name,
        assembler_args={
            "truncate_ratio_3p": truncate_ratio_3p,
            "truncate_ratio_5p": truncate_ratio_5p
        },
        adapter_args={
            "llrg_executable_path": llrg_executable_path,
            "hmm_model": hmm_model,
            "ccs_num_threads": 1,
            "samtools_executable_path": samtools_path,
            "ccs_executable_path": ccs_path,
            "ccs_pass": ccs_pass,
            "hmm_method": hmm_method,
            "strategy": strategy,
            "other_args": other_args
        },
        adapter_class=pbsim3.Pbsim3Adapter,
        is_pair_end=False
    )


def simulate_fast_mode(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        hmm_model: str,
        llrg_executable_path: str,
        depth: DepthType,
        jobs: int,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float,
        hmm_method: str,
        samtools_path: str,
        ccs_path: str,
        ccs_pass: int,
        strategy: str,
        simulator_name: Optional[str],
        other_args: List[str]
):
    ...


def main(args: List[str]):
    args, other_args = _parse_args(args)
    if args.ccs_pass > 1:
        samtools_path = enhanced_which(args.samtools_path)
        ccs_path = enhanced_which(args.ccs_path)
    else:
        samtools_path = ""
        ccs_path = ""
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        hmm_model=args.hmm_model,
        llrg_executable_path=enhanced_which(args.llrg_executable_path),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
        hmm_method=args.hmm_method,
        samtools_path=samtools_path,
        ccs_path=ccs_path, ccs_pass=args.ccs_pass,
        simulator_name=args.simulator_name,
        other_args=other_args,
        strategy=args.strategy
    )
