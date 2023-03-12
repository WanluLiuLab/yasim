import argparse
from typing import List, Tuple, Optional

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import patch_frontend_parser_public, enhanced_which
from yasim.llrg_adapter import art
from yasim.main import abstract_simulate

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--exename', required=False, help="Executable name or absolute path of art_illumina",
                        nargs='?',
                        type=str, action='store', default="art_illumina")
    parser.add_argument(
        "--sequencer",
        required=False,
        help="Name of Illumina Sequencer to Simulate: " + ", ".join((
            f"{name} -- {art.AVAILABLE_ILLUMINA_ART_SEQUENCER[name][0]}"
            for name in art.AVAILABLE_ILLUMINA_ART_SEQUENCER.keys()
        )),
        nargs='?',
        choices=art.AVAILABLE_ILLUMINA_ART_SEQUENCER.keys(),
        type=str,
        action='store',
        default="HS25"
    )
    parser.add_argument(
        "--rlen",
        required=False,
        help="Read length. Sequencer -- Read Length Table: " + ", ".join((
            f"{v[0]} -- {v[1]}"
            for v in art.AVAILABLE_ILLUMINA_ART_SEQUENCER.values()
        )),
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    parser.add_argument(
        '--mflen_mean',
        required=False,
        help="[PE Only] The mean size of DNA/RNA fragments for paired-end simulations",
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    parser.add_argument(
        '--mflen_std',
        required=False,
        help="[PE Only] The standard deviation of DNA/RNA fragment size for paired-end simulations.",
        nargs='?',
        type=int,
        action='store',
        default=0
    )
    parser.add_argument(
        '--is_pair_end',
        required=False,
        help="Whether to use Pair End (PE) Simulation",
        action='store_true'
    )
    parser = patch_frontend_parser_public(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        exename: str,
        depth: DepthType,
        jobs: int,
        sequencer: str,
        rlen: int,
        mflen_mean: int,
        mflen_std: int,
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
            "exename": exename,
            "other_args": other_args,
            "sequencer": sequencer,
            "rlen": rlen,
            "mflen_mean": mflen_mean,
            "mflen_std": mflen_std,
            "is_pair_end": is_pair_end
        },
        adapter_class=art.ArtAdapter,
        is_pair_end=is_pair_end
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    if args.is_pair_end:
        if args.mflen_mean * args.mflen_std == 0:
            _lh.error("Please set mflen_mean and mflen_std in PE simulation")
            return
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        exename=enhanced_which(args.exename),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        simulator_name=args.simulator_name,
        other_args=other_args,
        sequencer=args.sequencer,
        rlen=args.rlen,
        mflen_mean=args.mflen_mean,
        mflen_std=args.mflen_std,
        is_pair_end=args.is_pair_end
    )
