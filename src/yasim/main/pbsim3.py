import argparse
import glob
import os.path
from typing import List, Tuple, Optional

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper import parallel_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType, read_depth
from yasim.helper.llrg import pair_depth_info_with_transcriptome_fasta_filename, assemble_single_end, \
    patch_frontend_parser, \
    enhanced_which, AssembleSingleEnd, generate_callback
from yasim.llrg_adapter import pbsim3_transcriptome as pbsim3

logger = get_logger(__name__)

QSHMM_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename.replace("QSHMM-", ""))[0])
    for filename in glob.glob(os.path.join(pbsim3.PBSIM3_DIST, "QSHMM-*.model"))
]

ERRHMM_POSSIBLE_MODELS = [
    os.path.basename(os.path.splitext(filename.replace("ERRHMM-", ""))[0])
    for filename in glob.glob(os.path.join(pbsim3.PBSIM3_DIST, "ERRHMM-*.model"))
]


def _parse_args(args: List[str]) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m',
        '--hmm_model',
        required=True,
        help="Basename of HMM file. "
             f"If you select errhmm in hmm_method, it would be {ERRHMM_POSSIBLE_MODELS}"
             f"If you select qshmm in hmm_method, it would be {QSHMM_POSSIBLE_MODELS}",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-e',
        '--exename',
        required=False,
        help="Executable name of pbsim3, may be pbsim3 or pbsim.",
        nargs='?',
        type=str,
        action='store',
        default="pbsim3"
    )
    parser.add_argument(
        "-M",
        "--hmm_method",
        required=True,
        help="Whether to simulate using quality score (as PBSIM2) or error profile (new)",
        nargs='?',
        type=str,
        action='store',
        choices=["errhmm", "qshmm"]
    )
    parser.add_argument(
        "--ccs_pass",
        required=False,
        help="CCS Multipass Settings. Use 1 for CLR and others for CCS.",
        nargs='?',
        type=int,
        action='store',
        default=1
    )
    parser.add_argument(
        '--ccs_path',
        required=False,
        help="Executable name of ccs or pbccs. Omitted if ccs_pass == 1.",
        nargs='?',
        type=str,
        action='store',
        default="ccs"
    )
    parser.add_argument(
        '--samtools_path',
        required=False,
        help="Executable name of samtools. Omitted if ccs_pass == 1.",
        nargs='?',
        type=str,
        action='store',
        default="samtools"
    )
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
        hmm_method: str,
        samtools_path: str,
        ccs_path: str,
        ccs_pass: int,
        simulator_name: Optional[str],
        other_args: List[str]
):
    if simulator_name is None:
        simulator_name = "_".join(("pbsim3", hmm_model, "ccs" if ccs_pass > 1 else "clr"))
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = parallel_helper.ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(pair_depth_info_with_transcriptome_fasta_filename(transcriptome_fasta_dir, depth))
    assembler = AssembleSingleEnd(
        depth=depth,
        output_fastq_prefix=output_fastq_prefix,
        simulator_name=simulator_name,
        truncate_ratio_3p=truncate_ratio_3p,
        truncate_ratio_5p=truncate_ratio_5p,
        input_transcriptome_fasta_dir=transcriptome_fasta_dir
    )
    assembler.start()
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        if transcript_depth == 0:
            continue
        sim_thread = pbsim3.Pbsim3Adapter(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            hmm_model=hmm_model,
            exename=exename,
            other_args=other_args,
            ccs_num_threads=1,
            samtools_path=samtools_path,
            ccs_path=ccs_path,
            ccs_pass=ccs_pass,
            hmm_method=hmm_method
        )
        simulating_pool.append(sim_thread, callback=generate_callback(assembler, transcript_id))
    simulating_pool.start()
    simulating_pool.join()
    assembler.terminate()
    assembler.join()


def main(args: List[str]):
    args, other_args = _parse_args(args)
    if args.ccs_pass > 1:
        samtools_path = enhanced_which(args.samtools_path)
        ccs_path = enhanced_which(args.ccs_path)
    else:
        samtools_path = ""
        ccs_path = ""
    simulate(
        transcriptome_fasta_dir=args.fastas,
        output_fastq_prefix=args.out,
        hmm_model=args.hmm_model,
        exename=enhanced_which(args.exename),
        depth=read_depth(args.depth),
        jobs=args.jobs,
        truncate_ratio_3p=args.truncate_ratio_3p,
        truncate_ratio_5p=args.truncate_ratio_5p,
        hmm_method=args.hmm_method,
        samtools_path=samtools_path,
        ccs_path=ccs_path, ccs_pass=args.ccs_pass,
        simulator_name=args.simulator_name,
        other_args=other_args
    )
