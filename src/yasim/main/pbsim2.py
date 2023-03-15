import argparse
import glob
import os.path
from typing import List, Tuple, Optional

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_tgs, enhanced_which
from yasim.llrg_adapter import pbsim2
from yasim.main import abstract_simulate

logger = get_logger(__name__)

ALL_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename)[0])
    for filename in glob.glob(os.path.join(pbsim2.PBSIM2_DIST, "*.model"))
]


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--hmm_model', required=True, help="Basename of HMM file", nargs='?',
                        type=str, action='store', choices=ALL_POSSIBLE_MODELS)
    parser.add_argument('-e', '--exename', required=False,
                        help="Executable name of pbsim2, may be pbsim2 or pbsim.", nargs='?',
                        type=str, action='store', default="pbsim2")
    parser = patch_frontend_parser_tgs(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        exename: str,
        depth: DepthType,
        hmm_model: str,
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
        simulator_name="_".join(("pbsim2", hmm_model, "clr")) if simulator_name is None else simulator_name,
        assembler_args={
            "truncate_ratio_3p": truncate_ratio_3p,
            "truncate_ratio_5p": truncate_ratio_5p
        },
        adapter_args={
            "hmm_model": hmm_model,
            "exename": exename,
            "other_args": other_args
        },
        adapter_class=pbsim2.Pbsim2Adapter,
        is_pair_end=False
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        hmm_model=args.hmm_model,
        exename=enhanced_which(args.exename),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
        simulator_name=args.simulator_name,
        other_args=other_args
    )
