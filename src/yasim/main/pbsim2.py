import argparse
import glob
import os.path
from typing import List, Tuple, Optional

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper import parallel_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import pair_depth_info_with_transcriptome_fasta_filename, assemble_single_end, \
    patch_frontend_parser
from yasim.llrg_adapter import pbsim2

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
    parser = patch_frontend_parser(parser)
    return parser.parse_known_args(args)


def simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        hmm_model: str,
        exename: str,
        depth: DepthType,
        jobs: int,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float,
        simulator_name: Optional[str],
        other_args: List[str]
):
    if simulator_name is None:
        simulator_name = "_".join(("pbsim2", hmm_model, "CLR"))
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = parallel_helper.ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(pair_depth_info_with_transcriptome_fasta_filename(transcriptome_fasta_dir, depth))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        if transcript_depth == 0:
            continue
        sim_thread = pbsim2.Pbsim2Adapter(
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
        simulator_name=simulator_name,
        truncate_ratio_3p=truncate_ratio_3p,
        truncate_ratio_5p=truncate_ratio_5p,
        input_transcriptome_fasta_dir=transcriptome_fasta_dir
    )


def main(args: List[str]):
    args, other_args = _parse_args(args)
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        hmm_model=args.hmm_model,
        exename=args.exename,
        depth=read_depth(args.depth),
        jobs=args.jobs,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
        simulator_name=args.simulator_name,
        other_args=other_args
    )
