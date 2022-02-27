import argparse
import glob
import os.path
from typing import List

import commonutils.stdlib_helper.parallel_helper
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger
from yasim.main._helper import get_depth_from_intermediate_fasta, assemble_single_end
from yasim.simulator import badread

logger = get_logger(__name__)

__version__ = 0.1

ALL_POSSIBLE_MODELS = ["nanopore2018", "nanopore2020", "pacbio2016"]


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-m', '--model_name', required=True, help="Badread model name", nargs='?',
                        type=str, action='store', choices=ALL_POSSIBLE_MODELS)
    parser.add_argument('-e', '--badread_exename', required=False,
                        help="Executable name or absolute path of badread.", nargs='?',
                        type=str, action='store', default=None)
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))

    return parser.parse_args(args)


def simulate(
        intermediate_fasta_dir: str,
        output_fastq_prefix: str,
        model_name: str,
        badread_exename: str
):
    output_fastq_dir = output_fastq_prefix + ".d"
    simulating_pool = commonutils.stdlib_helper.parallel_helper.ParallelJobQueue(
        pool_name="Simulating jobs"
    )
    depth_info = list(get_depth_from_intermediate_fasta(intermediate_fasta_dir))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        sim_thread = badread.SimulatorBadread(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            model_name=model_name,
            badread_exename=badread_exename
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    assemble_single_end(depth_info, output_fastq_prefix, simulator_name=f"badread_{model_name}")


def main(args: List[str]):
    args = _parse_args(args)
    simulate(args.fastas, args.out, args.model_name, args.badread_exename)
