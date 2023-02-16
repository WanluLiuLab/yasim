import argparse
import multiprocessing
import os
import shutil
import stat
import threading
import time
from typing import Iterable, Tuple, List

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.parser.fastq import FastqWriter, FastqIterator
from labw_utils.bioutils.record.fastq import FastqRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io import file_system
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.depth import DepthType

DepthInfoType = Iterable[Tuple[int, str, str]]
"""
Depth information used by LLRG frontend interfaces.

They are: [depth, transcript_id, filename]
"""

_lh = get_logger(__name__)


def pair_depth_info_with_transcriptome_fasta_filename(
        input_transcriptome_fasta_dir: str,
        depth: DepthType
) -> DepthInfoType:
    """
    Glob and parse a filename line ``base_dir/1/transcript_id.fasta``.
    """
    for transcript_id, transcript_depth in depth.items():
        filename = os.path.join(input_transcriptome_fasta_dir, transcript_id + ".fa")
        yield transcript_depth, transcript_id, filename


def remark_fastq_single_end(
        input_filename: str,
        writer: FastqWriter,
        transcript_id: str,
        transcript_depth: float,
        simulator_name: str,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float
) -> Tuple[int, int]:
    """
    Re-mark all seq_id in FASTQ files.
    """
    num_of_reads = 0
    num_of_bases = 0
    for fastq_record in FastqIterator(input_filename, show_tqdm=False):
        sequence, quality = fastq_record.sequence, fastq_record.quality
        seq_len = len(sequence)
        truncate_len_3p = int(seq_len * truncate_ratio_3p)
        truncate_len_5p = int(seq_len * truncate_ratio_5p)
        sequence = sequence[truncate_len_3p:seq_len - truncate_len_5p]
        quality = quality[truncate_len_3p:seq_len - truncate_len_5p]
        new_record = FastqRecord(
            seq_id=f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}",
            sequence=sequence,
            quality=quality
        )
        writer.write(new_record)
        num_of_reads += 1
        num_of_bases += len(quality)
    return num_of_reads, num_of_bases


def remark_fastq_pair_end(
        input_filename_1: str,
        input_filename_2: str,
        writer1: FastqWriter,
        writer2: FastqWriter,
        transcript_id: str,
        transcript_depth: float,
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
        new_fastq_record_1 = FastqRecord(
            seq_id=f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}/1",
            sequence=fastq_record_1.sequence,
            quality=fastq_record_1.quality
        )
        new_fastq_record_2 = FastqRecord(
            seq_id=f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}/2",
            sequence=fastq_record_2.sequence,
            quality=fastq_record_2.quality
        )
        writer1.write(new_fastq_record_1)
        writer2.write(new_fastq_record_2)
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


class AssembleSingleEnd(threading.Thread):
    _transcript_ids_pending: List[str]
    _should_stop: bool

    def __init__(
            self,
            depth: DepthType,
            output_fastq_prefix: str,
            simulator_name: str,
            truncate_ratio_3p: float,
            truncate_ratio_5p: float,
            input_transcriptome_fasta_dir: str
    ):
        super().__init__()
        self._transcript_ids_pending = []
        self._should_stop = False
        self._depth = depth
        self._output_fastq_prefix = output_fastq_prefix
        self._simulator_name = simulator_name
        self._truncate_ratio_3p = truncate_ratio_3p
        self._truncate_ratio_5p = truncate_ratio_5p
        self._input_transcriptome_fasta_dir = input_transcriptome_fasta_dir

    def run(self):
        output_fastq_dir = self._output_fastq_prefix + ".d"
        with FastqWriter(self._output_fastq_prefix + ".fq") as writer, \
                get_writer(self._output_fastq_prefix + ".fq.stats") as stats_writer:
            stats_writer.write("\t".join((
                "TRANSCRIPT_ID",
                "INPUT_DEPTH",
                "SIMULATED_N_OF_READS",
                "SIMULATED_N_OF_BASES",
                "TRANSCRIBED_LENGTH",
                "SIMULATED_DEPTH"
            )) + "\n")

            while not self._should_stop:
                while not len(self._transcript_ids_pending) == 0:
                    transcript_id = self._transcript_ids_pending.pop(0)
                    this_fasta_path = os.path.join(self._input_transcriptome_fasta_dir, transcript_id + ".fa")
                    this_fastq_basename = os.path.join(output_fastq_dir, transcript_id)
                    this_fastq_path = this_fastq_basename + ".fq"
                    if not file_system.file_exists(this_fastq_path) or not file_system.file_exists(this_fasta_path):
                        _lh.error(f"Skipped error transcript: %s", transcript_id)
                        continue
                    try:
                        transcribed_length = FastaViewFactory(
                            filename=this_fasta_path,
                            read_into_memory=True,
                            show_tqdm=False
                        ).get_chr_length(transcript_id)
                    except KeyError:
                        _lh.error(f"Skipped error transcript: %s", transcript_id)
                        continue
                    transcript_depth = self._depth[transcript_id]
                    num_of_reads, num_of_bases = remark_fastq_single_end(
                        input_filename=this_fastq_path,
                        writer=writer,
                        transcript_id=transcript_id,
                        transcript_depth=transcript_depth,
                        simulator_name=self._simulator_name,
                        truncate_ratio_3p=self._truncate_ratio_3p,
                        truncate_ratio_5p=self._truncate_ratio_5p
                    )
                    stats_writer.write("\t".join((
                        transcript_id,  # "TRANSCRIPT_ID",
                        str(transcript_depth),  # "INPUT_DEPTH",
                        str(num_of_reads),  # "SIMULATED_N_OF_READS",
                        str(num_of_bases),  # "SIMULATED_N_OF_BASES",
                        str(transcribed_length),  # "TRANSCRIBED_LENGTH",
                        str(num_of_bases / transcribed_length)  # "SIMULATED_DEPTH"
                    )) + "\n")
                time.sleep(0.01)

    def add_transcript_id(self, transcript_id: str):
        self._transcript_ids_pending.append(transcript_id)

    def terminate(self):
        self._should_stop = True


def assemble_single_end(
        depth: DepthType,
        output_fastq_prefix: str,
        simulator_name: str,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float,
        input_transcriptome_fasta_dir: str
):
    """
    Assemble single_end reads into one.
    """
    output_fastq_dir = output_fastq_prefix + ".d"
    with FastqWriter(output_fastq_prefix + ".fq") as writer, \
            get_writer(output_fastq_prefix + ".fq.stats") as stats_writer:
        stats_writer.write("\t".join((
            "TRANSCRIPT_ID",
            "INPUT_DEPTH",
            "SIMULATED_N_OF_READS",
            "SIMULATED_N_OF_BASES",
            "TRANSCRIBED_LENGTH",
            "SIMULATED_DEPTH"
        )) + "\n")
        for transcript_id, transcript_depth in tqdm(iterable=depth.items(), desc="Merging..."):
            this_fasta_path = os.path.join(input_transcriptome_fasta_dir, transcript_id + ".fa")
            this_fastq_basename = os.path.join(output_fastq_dir, transcript_id)
            this_fastq_path = this_fastq_basename + ".fq"
            if not file_system.file_exists(this_fastq_path) or not file_system.file_exists(this_fasta_path):
                _lh.error(f"Skipped error transcript: %s", transcript_id)
                continue
            try:
                transcribed_length = FastaViewFactory(
                    filename=this_fasta_path,
                    read_into_memory=True,
                    show_tqdm=False
                ).get_chr_length(transcript_id)
            except KeyError:
                _lh.error(f"Skipped error transcript: %s", transcript_id)
                continue

            num_of_reads, num_of_bases = remark_fastq_single_end(
                input_filename=this_fastq_path,
                writer=writer,
                transcript_id=transcript_id,
                transcript_depth=transcript_depth,
                simulator_name=simulator_name,
                truncate_ratio_3p=truncate_ratio_3p,
                truncate_ratio_5p=truncate_ratio_5p
            )
            stats_writer.write("\t".join((
                transcript_id,  # "TRANSCRIPT_ID",
                str(transcript_depth),  # "INPUT_DEPTH",
                str(num_of_reads),  # "SIMULATED_N_OF_READS",
                str(num_of_bases),  # "SIMULATED_N_OF_BASES",
                str(transcribed_length),  # "TRANSCRIBED_LENGTH",
                str(num_of_bases / transcribed_length)  # "SIMULATED_DEPTH"
            )) + "\n")


def patch_frontend_parser(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument('-F', '--fastas', required=True,
                        help="Directory of transcribed DGE FASTAs from `transcribe` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-d', '--depth', required=True, help="Depth generated by `dge` step", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output transcript prefix", nargs='?',
                        type=str, action='store')
    parser.add_argument('-j', '--jobs', required=False,
                        help="Number of threads", nargs='?',
                        type=int, action='store', default=multiprocessing.cpu_count())
    parser.add_argument('--truncate_ratio_3p', required=False,
                        help="Ratio of 3 prime truncation", nargs='?',
                        type=float, action='store', default=0.0)
    parser.add_argument('--truncate_ratio_5p', required=False,
                        help="Ratio of 5 prime truncation", nargs='?',
                        type=float, action='store', default=0.0)
    parser.add_argument('--simulator_name', required=False,
                        help="Custom simulator name. Used in FASTQ tags", nargs='?',
                        type=str, action='store', default=None)
    return parser


def enhanced_which(path_or_filename: str) -> str:
    if os.path.exists(path_or_filename):
        stat_result = os.stat(path_or_filename)
        if not stat.S_IXUSR & stat_result.st_mode:
            raise FileNotFoundError(f"File {path_or_filename} is not executable!")
    else:
        resolved_path = shutil.which(path_or_filename)
        if resolved_path is None:
            raise FileNotFoundError(f"File {path_or_filename} not found!")
        else:
            path_or_filename = resolved_path
    return path_or_filename
