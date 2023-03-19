import argparse
from typing import List, Tuple, Optional

from yasim._main import abstract_simulate
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_tgs, enhanced_which, patch_frontend_parser_public
from yasim.llrg_adapter import pbsim


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser = patch_frontend_parser_public(
        parser,
        llrg_name="pbsim",
        default_llrg_executable_name="pbsim"
    )
    parser = patch_frontend_parser_tgs(parser)
    parser = pbsim.patch_frontend_parser(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        llrg_executable_path: str,
        depth: DepthType,
        is_ccs: bool,
        jobs: int,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float,
        simulator_name: Optional[str],
        other_args: List[str]
):
    abstract_simulate(
        transcriptome_fasta_dir=transcriptome_fasta_dir,
        output_fastq_prefix=output_fastq_prefix,
        depth=depth,
        jobs=jobs,
        simulator_name="_".join(
            ("pbsim", "RS", "ccs" if is_ccs else "clr")) if simulator_name is None else simulator_name,
        assembler_args={
            "truncate_ratio_3p": truncate_ratio_3p,
            "truncate_ratio_5p": truncate_ratio_5p
        },
        adapter_args={
            "is_ccs": is_ccs,
            "llrg_executable_path": llrg_executable_path,
            "other_args": other_args
        },
        adapter_class=pbsim.PbsimAdapter,
        is_pair_end=False
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    depth = read_depth(args.depth)
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        llrg_executable_path=enhanced_which(args.llrg_executable_path),
        depth=depth,
        is_ccs=args.ccs,
        jobs=args.jobs,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
        simulator_name=args.simulator_name,
        other_args=other_args
    )