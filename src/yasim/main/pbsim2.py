import argparse
import glob
import os.path
from typing import List

import commonutils.parallel_helper
from commonutils.logger import get_logger
from commonutils.tqdm_importer import tqdm
from yasim.main._helper import get_depth_from_intermediate_fasta
from yasim.simulator import pbsim2

logger = get_logger(__name__)

__version__ = 0.1

ALL_POSSIBLE_MODELS = [os.path.basename(os.path.splitext(filename)[0]) for filename in
                       glob.glob(os.path.join(pbsim2.FILE_DIR, "pbsim2_dist", f"*.model"))]


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-m', '--hmm_model', required=True, help="Basename of HMM file", nargs='?',
                        type=str, action='store', choices=ALL_POSSIBLE_MODELS)
    parser.add_argument('-e', '--pbsim2_exename', required=False,
                        help="Executable name of pbsim2, may be pbsim2 or pbsim.", nargs='?',
                        type=str, action='store', default="pbsim2")
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))

    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        hmm_model: str,
        pbsim2_exename: str
):
    output_fastq_dir = output_fastq_prefix + ".d"
    simulating_pool = commonutils.parallel_helper.ParallelJobQueue(
        pool_name="Simulating jobs"
    )
    transcript_depths = get_depth_from_intermediate_fasta(intermediate_fasta_dir)
    for transcript_depth in tqdm(iterable=transcript_depths, desc="Submitting jobs..."):
        sim_thread = pbsim2.SimulatePbsim2(
            input_fasta=os.path.join(intermediate_fasta_dir, f"{transcript_depth}.fa"),
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_depth),
            depth=transcript_depth,
            hmm_model=hmm_model,
            pbsim2_exename=pbsim2_exename
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    with open(output_fastq_prefix + ".fq", "wb") as writer:
        for file in tqdm(iterable=list(glob.glob(os.path.join(output_fastq_dir, "*.fq"))),
                         desc="Assembling read"):
            writer.write(open(file, "rb").read())


def main(args: List[str]):
    args = _parse_args(args)
    simulate(args.fastas, args.out, args.hmm_model, args.pbsim2_exename)
