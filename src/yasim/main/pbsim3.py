import argparse
import glob
import multiprocessing
import os.path
from typing import List, Tuple

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper import parallel_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import get_depth_from_intermediate_fasta, assemble_single_end
from yasim.llrg_adapter import pbsim3_transcriptome as pbsim3

logger = get_logger(__name__)

ALL_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename.replace("QSHMM-", ""))[0])
    for filename in glob.glob(os.path.join(pbsim3.PBSIM3_DIST, "QSHMM-*.model"))
]


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--depth', required=True, help="Depth generated by `dge` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-m', '--hmm_model', required=True, help="Basename of HMM file", nargs='?',
                        type=str, action='store', choices=ALL_POSSIBLE_MODELS)
    parser.add_argument('-e', '--exename', required=False,
                        help="Executable name of pbsim3, may be pbsim3 or pbsim.", nargs='?',
                        type=str, action='store', default="pbsim3")
    parser.add_argument('-j', '--jobs', required=False,
                        help="Number of threads", nargs='?',
                        type=int, action='store', default=multiprocessing.cpu_count())
    return parser.parse_known_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        hmm_model: str,
        exename: str,
        depth: DepthType,
        jobs: int,
        other_args:List[str]
):
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = parallel_helper.ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(get_depth_from_intermediate_fasta(intermediate_fasta_dir, depth))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        sim_thread = pbsim3.Pbsim3Adapter(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            hmm_model=hmm_model,
            exename=exename,
            other_args=other_args
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    assemble_single_end(
        depth=depth,
        output_fastq_prefix=output_fastq_prefix,
        simulator_name=f"pbsim3_{hmm_model}"
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    depth = read_depth(args.depth)
    simulate(
        intermediate_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        hmm_model=args.hmm_model,
        exename=args.exename,
        depth=depth,
        jobs=args.jobs,
        other_args=other_args
    )
