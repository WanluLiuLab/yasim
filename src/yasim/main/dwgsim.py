import argparse
import os.path
import warnings
from typing import List, Optional

import commonutils.stdlib_helper.parallel_helper
import yasim.simulator.dwgsim
from commonutils.importer.tqdm_importer import tqdm
from commonutils.stdlib_helper.logger_helper import get_logger
from yasim import simulator
from yasim.main._helper import get_depth_from_intermediate_fasta, assemble_pair_end

warnings.warn("NGS simulator is not used, so not updated.", DeprecationWarning, stacklevel=2)

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
    simulating_pool = commonutils.stdlib_helper.parallel_helper.ParallelJobQueue(
        pool_name="Simulating jobs"
    )
    depth_info = list(get_depth_from_intermediate_fasta(intermediate_fasta_dir))
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        sim_thread = simulator.dwgsim.SimulatorDwgsim(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            dwgsim_exename=dwgsim_exename
        )
        simulating_pool.append(sim_thread)
    simulating_pool.start()
    simulating_pool.join()
    assemble_pair_end(depth_info, output_fastq_prefix, simulator_name="dwgsim")



def main(args: List[str]):
    args = _parse_args(args)
    simulate(args.fastas, args.out, args.dwgsim_exename)
