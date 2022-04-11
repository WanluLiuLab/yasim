import argparse
import multiprocessing
import os.path
from typing import List

from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger
from commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import get_depth_from_intermediate_fasta, assemble_single_end
from yasim.llrg_adapter import badread

logger = get_logger(__name__)

ALL_POSSIBLE_MODELS = ("nanopore2018", "nanopore2020", "pacbio2016", "verybad", "verynice")
"""All possible badread model names"""


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--depth', required=True, help="Depth generated by `dge` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-m', '--model_name', required=True, help="Badread model name", nargs='?',
                        type=str, action='store', choices=ALL_POSSIBLE_MODELS)
    parser.add_argument('-e', '--exename', required=False,
                        help="Executable name or absolute path of badread.", nargs='?',
                        type=str, action='store', default=None)
    parser.add_argument('-j', '--jobs', required=False,
                        help="Number of threads", nargs='?',
                        type=int, action='store', default=multiprocessing.cpu_count())
    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        model_name: str,
        exename: str,
        depth: DepthType,
        jobs: int
):
    output_fastq_dir = output_fastq_prefix + ".d"
    simulating_pool = ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(get_depth_from_intermediate_fasta(intermediate_fasta_dir, depth))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        sim_thread = badread.BadreadAdapter(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            model_name=model_name,
            exename=exename
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    assemble_single_end(depth, output_fastq_prefix, simulator_name=f"badread_{model_name}")


def main(args: List[str]):
    args = _parse_args(args)
    depth = read_depth(args.depth)
    simulate(
        intermediate_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        model_name=args.model_name,
        exename=args.exename,
        depth=depth,
        jobs=args.jobs
    )
