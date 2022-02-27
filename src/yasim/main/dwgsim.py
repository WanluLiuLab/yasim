import argparse
import os.path
from typing import List, Optional

import commonutils.stdlib_helper.parallel_helper
import yasim.simulator.dwgsim
from bioutils.io.fastq import FastqWriter
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer
from commonutils.stdlib_helper.logger_helper import get_logger
from yasim import simulator
from yasim.main._helper import remark_fastq_pair_end, get_depth_from_intermediate_fasta

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
    with FastqWriter(output_fastq_prefix + "_1.fq") as writer1, \
            FastqWriter(output_fastq_prefix + "_2.fq") as writer2, \
            get_writer(output_fastq_prefix + ".fq.stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "THEORETICAL_DEPTH",
            "ACTUAL_N_OF_READS",
        ))+"\n")
        for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Merging..."):
            this_fastq_basename = os.path.join(output_fastq_dir, transcript_id)
            num_of_reads = remark_fastq_pair_end(
                input_filename_1=this_fastq_basename+"_1.fq",
                input_filename_2=this_fastq_basename + "_2.fq",
                writer1=writer1,
                writer2=writer2,
                transcript_id=transcript_id,
                transcript_depth=transcript_depth,
                simulator_name="dwgsim"
            )
            stats_writer.write("\t".join((
                transcript_id,
                str(transcript_depth),
                str(num_of_reads)
            ))+"\n")


def main(args: List[str]):
    args = _parse_args(args)
    simulate(args.fastas, args.out, args.dwgsim_exename)
