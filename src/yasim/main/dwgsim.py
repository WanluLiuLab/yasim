import argparse
from typing import List, Tuple, Optional

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_public, enhanced_which
from yasim.llrg_adapter import dwgsim
from yasim.main import abstract_simulate

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--exename', required=False, help="Executable name or absolute path of dwgsim",
                        nargs='?',
                        type=str, action='store', default="dwgsim")
    parser = patch_frontend_parser_public(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        exename: str,
        depth: DepthType,
        jobs: int,
        simulator_name: Optional[str],
        other_args: List[str]
):
    abstract_simulate(
        transcriptome_fasta_dir=transcriptome_fasta_dir,
        output_fastq_prefix=output_fastq_prefix,
        depth=depth,
        jobs=jobs,
        simulator_name="dwgsim" if simulator_name is None else simulator_name,
        assembler_args={},
        adapter_args={
            "exename": exename,
            "other_args": other_args
        },
        adapter_class=dwgsim.DwgsimAdapter,
        is_pair_end=True
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        exename=enhanced_which(args.exename),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        simulator_name=args.simulator_name,
        other_args=other_args
    )
