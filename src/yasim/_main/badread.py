import argparse
from typing import List, Tuple, Optional

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim._main import abstract_simulate
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_tgs, enhanced_which, patch_frontend_parser_public
from yasim.llrg_adapter import badread

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser = patch_frontend_parser_tgs(parser)
    parser = patch_frontend_parser_public(parser, "badread", "badread")
    parser = badread.patch_frontend_parser(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        model_name: str,
        llrg_executable_path: str,
        depth: DepthType,
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
        simulator_name="_".join(("badread", model_name, "cDNA")) if simulator_name is None else simulator_name,
        assembler_args={
            "truncate_ratio_3p": truncate_ratio_3p,
            "truncate_ratio_5p": truncate_ratio_5p
        },
        adapter_args={
            "model_name": model_name,
            "llrg_executable_path": llrg_executable_path,
            "other_args": other_args
        },
        adapter_class=badread.BadReadAdapter,
        is_pair_end=False
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        model_name=args.model_name,
        llrg_executable_path=enhanced_which(args.llrg_executable_path),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
        simulator_name=args.simulator_name,
        other_args=other_args
    )
