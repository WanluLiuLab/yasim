import os
from typing import Iterable, Tuple

from bioutils.io.fastq import FastqWriter, FastqIterator
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io import file_system
from commonutils.io.safe_io import get_writer
from yasim.helper.depth import DepthType

DepthInfoType = Iterable[Tuple[int, str, str]]
"""
Depth information used by LLRG frontend interfaces.

They are: [depth, transcript_id, filename]
"""


def get_depth_from_intermediate_fasta(
        intermediate_fasta_dir: str,
        depth: DepthType
) -> DepthInfoType:
    """
    Glob and parse a filename line base_dir/1/transcript_id.fasta.
    """
    for transcript_id, transcript_depth in depth.items():
        filename = os.path.join(intermediate_fasta_dir, transcript_id + ".fa")
        yield transcript_depth, transcript_id, filename


def remark_fastq_single_end(
        input_filename: str,
        writer: FastqWriter,
        transcript_id: str,
        transcript_depth: int,
        simulator_name: str
) -> int:
    """
    Re-mark all seq_id in FASTQ files.
    """
    num_of_reads = 0
    for fastq_record in FastqIterator(input_filename, show_tqdm=False):
        fastq_record.seq_id = f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}"
        writer.write(fastq_record)
        num_of_reads += 1
    return num_of_reads


def remark_fastq_pair_end(
        input_filename_1: str,
        input_filename_2: str,
        writer1: FastqWriter,
        writer2: FastqWriter,
        transcript_id: str,
        transcript_depth: int,
        simulator_name: str
) -> int:
    """
    Re-mark all seq_id in FASTQ files, return number of reads
    """
    num_of_reads = 0
    for fastq_record_1, fastq_record_2 in zip(
            FastqIterator(input_filename_1, show_tqdm=False),
            FastqIterator(input_filename_2, show_tqdm=False)
    ):
        fastq_record_1.seq_id = f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}/1"
        fastq_record_2.seq_id = f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}/2"
        writer1.write(fastq_record_1)
        writer2.write(fastq_record_2)
        num_of_reads += 1
    return num_of_reads


def assemble_pair_end(
        depth: DepthType,
        output_fastq_prefix: str,
        simulator_name: str
):
    """
    Assemble pair-end reads into one.
    """
    output_fastq_dir = output_fastq_prefix + ".d"
    with FastqWriter(output_fastq_prefix + "_1.fq") as writer1, \
            FastqWriter(output_fastq_prefix + "_2.fq") as writer2, \
            get_writer(output_fastq_prefix + ".fq.stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "INPUT_DEPTH",
            "SIMULATED_N_OF_READS",
        )) + "\n")
        for transcript_id, transcript_depth in tqdm(iterable=depth.items(), desc="Merging..."):
            this_fastq_basename = os.path.join(output_fastq_dir, transcript_id)
            num_of_reads = remark_fastq_pair_end(
                input_filename_1=this_fastq_basename + "_1.fq",
                input_filename_2=this_fastq_basename + "_2.fq",
                writer1=writer1,
                writer2=writer2,
                transcript_id=transcript_id,
                transcript_depth=transcript_depth,
                simulator_name=simulator_name
            )
            stats_writer.write("\t".join((
                transcript_id,
                str(transcript_depth),
                str(num_of_reads)
            )) + "\n")


def assemble_single_end(
        depth: DepthType,
        output_fastq_prefix: str,
        simulator_name: str
):
    """
    Assemble single_end reads into one.
    """
    output_fastq_dir = output_fastq_prefix + ".d"
    with FastqWriter(output_fastq_prefix + ".fq") as writer, get_writer(
            output_fastq_prefix + ".fq.stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "INPUT_DEPTH",
            "SIMULATED_N_OF_READS",
        )) + "\n")
        for transcript_id, transcript_depth in tqdm(iterable=depth.items(), desc="Merging..."):
            this_fastq_basename = os.path.join(output_fastq_dir, transcript_id)
            if not file_system.file_exists(this_fastq_basename + ".fq"):
                continue
            num_of_reads = remark_fastq_single_end(
                input_filename=this_fastq_basename + ".fq",
                writer=writer,
                transcript_id=transcript_id,
                transcript_depth=transcript_depth,
                simulator_name=simulator_name
            )
            stats_writer.write("\t".join((
                transcript_id,
                str(transcript_depth),
                str(num_of_reads)
            )) + "\n")
