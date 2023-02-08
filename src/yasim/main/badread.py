import argparse
import os.path
from typing import List, Tuple

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor

from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import get_depth_from_intermediate_fasta, assemble_single_end, patch_frontend_parser
from yasim.llrg_adapter import badread

logger = get_logger(__name__)

ALL_POSSIBLE_MODELS = ("nanopore2018", "nanopore2020", "pacbio2016", "verybad", "verynice")
"""All possible badread model names"""


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model_name', required=True, help="Badread model name", nargs='?',
                        type=str, action='store', choices=ALL_POSSIBLE_MODELS)
    parser.add_argument('-e', '--exename', required=False,
                        help="Executable name or absolute path of badread.", nargs='?',
                        type=str, action='store', default="badread")
    parser = patch_frontend_parser(parser)
    return parser.parse_known_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        model_name: str,
        exename: str,
        depth: DepthType,
        jobs: int,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float,
        other_args: List[str]
):
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(get_depth_from_intermediate_fasta(intermediate_fasta_dir, depth))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        if transcript_depth == 0:
            continue
        sim_thread = badread.BadReadAdapter(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            model_name=model_name,
            exename=exename,
            other_args=other_args
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    assemble_single_end(
        depth=depth,
        output_fastq_prefix=output_fastq_prefix,
        simulator_name=f"badread_{model_name}",
        truncate_ratio_3p=truncate_ratio_3p,
        truncate_ratio_5p=truncate_ratio_5p
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    depth = read_depth(args.depth)
    simulate(
        intermediate_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        model_name=args.model_name,
        exename=args.exename,
        depth=depth,
        jobs=args.jobs,
        other_args=other_args,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p
    )
