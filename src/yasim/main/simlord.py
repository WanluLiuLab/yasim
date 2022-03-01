import argparse
import glob
import os.path
import warnings
from typing import List

import commonutils.stdlib_helper.parallel_helper
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger
from yasim.main._helper import get_depth_from_intermediate_fasta
from yasim.simulator import simlord

warnings.warn("SimLoRD simulator is not used, so not updated.", DeprecationWarning, stacklevel=2)


logger = get_logger(__name__)

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-e', '--simlord_exename', required=False,
                        help="Executable name or absolute path of simlord.", nargs='?',
                        type=str, action='store', default=None)
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))

    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        simlord_exename: str
):
    output_fastq_dir = output_fastq_prefix + ".d"
    simulating_pool = commonutils.stdlib_helper.parallel_helper.ParallelJobQueue(
        pool_name="Simulating jobs"
    )
    transcript_depths = get_depth_from_intermediate_fasta(intermediate_fasta_dir)
    for transcript_depth in tqdm(iterable=transcript_depths, desc="Submitting jobs..."):
        sim_thread = simlord.SimulatorSimlord(
            input_fasta=os.path.join(intermediate_fasta_dir, f"{transcript_depth}.fa"),
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_depth),
            depth=transcript_depth,
            simlord_exename=simlord_exename
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    with open(output_fastq_prefix + ".fq", "wb") as writer:
        for file in tqdm(iterable=list(glob.glob(os.path.join(output_fastq_dir, "*.fastq"))),
                         desc="Assembling read"):
            writer.write(open(file, "rb").read())


def main(args: List[str]):
    args = _parse_args(args)
    simulate(args.fastas, args.out, args.simlord_exename)
