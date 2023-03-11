import argparse
import os.path
from typing import List, Tuple, Optional

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper import parallel_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import generate_callback, pair_depth_info_with_transcriptome_fasta_filename, AssemblePairEnd, \
    patch_frontend_parser
from yasim.llrg_adapter import dwgsim

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--exename', required=False, help="Executable name or absolute path of dwgsim",
                        nargs='?',
                        type=str, action='store', default="dwgsim")
    parser = patch_frontend_parser(parser)
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
    if simulator_name is None:
        simulator_name = "dwgsim"
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = parallel_helper.ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(pair_depth_info_with_transcriptome_fasta_filename(transcriptome_fasta_dir, depth))
    assembler = AssemblePairEnd(
        depth=depth,
        output_fastq_prefix=output_fastq_prefix,
        simulator_name=simulator_name,
        input_transcriptome_fasta_dir=transcriptome_fasta_dir
    )
    assembler.start()
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        if transcript_depth == 0 or not os.path.exists(transcript_filename):
            continue
        sim_thread = dwgsim.DwgsimAdapter(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            exename=exename,
            other_args=other_args
        )
        simulating_pool.append(sim_thread, callback=generate_callback(assembler, transcript_id))
    simulating_pool.start()
    simulating_pool.join()
    assembler.terminate()
    assembler.join()


def main(args: List[str]):
    args, other_args = _parse_args(args)
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        exename=args.exename,
        depth=read_depth(args.depth),
        jobs=args.jobs,
        simulator_name=args.simulator_name,
        other_args=other_args
    )
