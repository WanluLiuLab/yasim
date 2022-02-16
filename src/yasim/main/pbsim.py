import argparse
import glob
import os.path
from typing import List, Optional

import commonutils.parallel_helper
from commonutils import ioctl
from commonutils.logger import get_logger
from commonutils.tqdm_importer import tqdm
from yasim.main._helper import get_depth_from_intermediate_fasta
from yasim.simulator import pbsim

logger = get_logger(__name__)

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-c', '--ccs', required=False, help="Simulate CCS instead of CLR", action='store_true')
    parser.add_argument('-e', '--pbsim_exename', required=False, help="Executable name or absolute path of pbsim",
                        nargs='?',
                        type=str, action='store', default=None)
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        is_ccs: bool = False,
        pbsim_exename: Optional[str] = None
):
    output_fastq_dir = output_fastq_prefix + ".d"
    ioctl.mkdir_p(output_fastq_dir)
    simulating_pool = commonutils.parallel_helper.ParallelJobQueue(
        pool_name="Simulating jobs"
    )
    transcript_depths = get_depth_from_intermediate_fasta(intermediate_fasta_dir)
    for transcript_depth in tqdm(iterable=transcript_depths, desc="Submitting jobs..."):
        sim_thread = pbsim.SimulatorPbsim(
            input_fasta=os.path.join(intermediate_fasta_dir, f"{transcript_depth}.fa"),
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_depth),
            depth=transcript_depth,
            pbsim_exename=pbsim_exename,
            is_ccs=is_ccs
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
    simulate(args.fastas, args.out, args.ccs, args.pbsim_exename)
