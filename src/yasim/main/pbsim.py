import argparse
import multiprocessing
import os.path
from typing import List

import commonutils.io.file_system
import commonutils.shell_utils
import commonutils.stdlib_helper.parallel_helper
from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import get_depth_from_intermediate_fasta, assemble_single_end
from yasim.llrg_adapter import pbsim

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--depth', required=True, help="Depth generated by `dge` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-c', '--ccs', required=False, help="Simulate CCS instead of CLR", action='store_true')
    parser.add_argument('-e', '--exename', required=False, help="Executable name or absolute path of pbsim",
                        nargs='?',
                        type=str, action='store', default="pbsim")
    parser.add_argument('-j', '--jobs', required=False,
                        help="Number of threads", nargs='?',
                        type=int, action='store', default=multiprocessing.cpu_count)
    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        exename:str,
        depth:DepthType,
        is_ccs: bool,
        jobs:int
):
    output_fastq_dir = output_fastq_prefix + ".d"
    shell_utils.mkdir_p(output_fastq_dir)
    simulating_pool = commonutils.stdlib_helper.parallel_helper.ParallelJobExecutor(
        pool_name="Simulating jobs"
    )
    depth_info = list(get_depth_from_intermediate_fasta(intermediate_fasta_dir, depth))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        sim_thread = pbsim.PbsimAdapter(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            exename=exename,
            is_ccs=is_ccs
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    simulator_name = "pbsim_"
    if is_ccs:
        simulator_name += "ccs"
    else:
        simulator_name += "clr"
    assemble_single_end(depth, output_fastq_prefix, simulator_name=simulator_name)


def main(args: List[str]):
    args = _parse_args(args)
    depth=read_depth(args.depth)
    simulate(
        intermediate_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        exename=args.exename,
        depth=depth,
        is_ccs=args.ccs,
        jobs=args.jobs
    )
