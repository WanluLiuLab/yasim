import argparse
import glob
import os.path
from typing import List, Optional

import commonutils.parallel_helper
import yasim.simulator.dwgsim
from commonutils.logger import get_logger
from commonutils.tqdm_importer import tqdm
from yasim import simulator
from yasim.main._helper import get_depth_from_intermediate_fasta

logger = get_logger(__name__)

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    parser.add_argument('-e', '--dwgsim_exename', required=False, help="Executable name or absolute path of dwgsim",
                        nargs='?',
                        type=str, action='store', default=None)
    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        dwgsim_exename: Optional[str] = None
):
    output_fastq_dir = output_fastq_prefix + ".d"
    simulating_pool = commonutils.parallel_helper.ParallelJobQueue(
        pool_name="Simulating jobs"
    )
    transcript_depths = get_depth_from_intermediate_fasta(intermediate_fasta_dir)
    for transcript_depth in tqdm(iterable=transcript_depths, desc="Submitting jobs..."):
        sim_thread = simulator.dwgsim.SimulatorDwgsim(
            input_fasta=os.path.join(intermediate_fasta_dir, f"{transcript_depth}.fa"),
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_depth),
            depth=transcript_depth,
            dwgsim_exename=dwgsim_exename
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    with open(output_fastq_prefix + "_1.fq", "wb") as writer:
        for file in tqdm(iterable=list(glob.glob(os.path.join(output_fastq_dir, "*_1.fq"))),
                         desc="Assembling read1"):
            writer.write(open(file, "rb").read())
    with open(output_fastq_prefix + "_2.fq", "wb") as writer:
        for file in tqdm(iterable=list(glob.glob(os.path.join(output_fastq_dir, "*_2.fq"))),
                         desc="Assembling read2"):
            writer.write(open(file, "rb").read())


def main(args: List[str]):
    args = _parse_args(args)
    simulate(args.fastas, args.out, args.dwgsim_exename)
