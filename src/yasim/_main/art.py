__all__ = (
    "main",
    "simulate"
)

import argparse
from typing import List, Tuple, Optional

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim._main import abstract_simulate
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_public, enhanced_which
from yasim.llrg_adapter import art

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser = patch_frontend_parser_public(
        parser,
        llrg_name="art",
        default_llrg_executable_name="art_illumina"
    )
    parser = art.patch_frontend_parser(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        llrg_executable_path: str,
        depth: DepthType,
        jobs: int,
        sequencer_name: str,
        read_length: int,
        pair_end_fragment_length_mean: int,
        pair_end_fragment_length_std: int,
        is_pair_end: bool,
        simulator_name: Optional[str],
        other_args: List[str]
):
    abstract_simulate(
        transcriptome_fasta_dir=transcriptome_fasta_dir,
        output_fastq_prefix=output_fastq_prefix,
        depth=depth,
        jobs=jobs,
        simulator_name="art_illumina" if simulator_name is None else simulator_name,
        assembler_args={},
        adapter_args={
            "llrg_executable_path": llrg_executable_path,
            "other_args": other_args,
            "sequencer_name": sequencer_name,
            "read_length": read_length,
            "pair_end_fragment_length_mean": pair_end_fragment_length_mean,
            "pair_end_fragment_length_std": pair_end_fragment_length_std,
            "is_pair_end": is_pair_end
        },
        adapter_class=art.ArtAdapter,
        is_pair_end=is_pair_end
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    if args.is_pair_end:
        if args.pair_end_fragment_length_mean * args.pair_end_fragment_length_std == 0:
            _lh.error("Please set pair_end_fragment_length_mean and pair_end_fragment_length_std in PE simulation")
            return
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        llrg_executable_path=enhanced_which(args.llrg_executable_path),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        simulator_name=args.simulator_name,
        other_args=other_args,
        sequencer_name=args.sequencer_name,
        read_length=args.read_length,
        pair_end_fragment_length_mean=args.pair_end_fragment_length_mean,
        pair_end_fragment_length_std=args.pair_end_fragment_length_std,
        is_pair_end=args.is_pair_end
    )
