import glob
import os
from typing import List

from bioutils.io.fastq import FastqIterator, FastqWriter


def get_depth_from_intermediate_fasta(intermediate_fasta_dir: str) -> List[str]:
    """
    For a filename line base_dir/1/transcript_id.fasta, return iterator of [depth, transcript_id, filename]
    """
    retl = []
    for filename in glob.glob(os.path.join(intermediate_fasta_dir, "*.fa")):
        retl.append(os.path.basename(os.path.splitext(filename)[0]))
    return retl

def remark_fastq_single_end(
        input_filename:str,
        output_filename: str,
        transcript_id:str,
        transcript_depth:int,
        simulator_name:str
):
    """
    Re-mark all seq_id in FASTQ files.
    """
    i = 0
    with FastqWriter(output_filename) as writer:
        for fastq_record in FastqIterator(input_filename):
            fastq_record.seq_id=f"{transcript_id}:{i}:{transcript_depth}:{simulator_name}"
            writer.write(fastq_record)


def remark_fastq_pair_end(
        input_filename_1:str,
        input_filename_2:str,
        output_filename_1: str,
        output_filename_2: str,
        transcript_id:str,
        transcript_depth:int,
        simulator_name:str
):
    """
    Re-mark all seq_id in FASTQ files.
    """
    i = 0
    with FastqWriter(output_filename_1) as writer1, FastqWriter(output_filename_2) as writer2:
        for fastq_record_1, fastq_record_2 in zip(
                FastqIterator(input_filename_1),
                FastqIterator(input_filename_2, show_tqdm=False)
        ) :
            fastq_record_1.seq_id=f"{transcript_id}:{i}:{transcript_depth}:{simulator_name}/1"
            fastq_record_2.seq_id = f"{transcript_id}:{i}:{transcript_depth}:{simulator_name}/2"
            writer1.write(fastq_record_1)
            writer2.write(fastq_record_2)
