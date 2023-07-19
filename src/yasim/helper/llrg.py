"""
llrg.py -- Helper functions in Low-Level Read Generators.

.. versionadded:: 3.1.5
"""

__all__ = (
    "patch_frontend_parser_tgs",
    "patch_frontend_parser_public",
    "patch_frontend_parser_bulk_rna_seq",
    "patch_frontend_parser_sc_rna_seq",
    "generate_callback",
    "BaseAssembler",
    "AssembleDumb",
    "AssemblePairEnd",
    "AssembleSingleEnd",
    "enhanced_which"
)

import argparse
import multiprocessing
import os
import shutil
import threading
import time
from abc import ABC, abstractmethod

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.parser.fastq import FastqWriter, FastqIterator
from labw_utils.bioutils.record.fastq import FastqRecord, MisFormattedFastqRecordError
from labw_utils.commonutils.lwio import file_system
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Tuple, List, Callable, Optional
from yasim.helper.depth_io import DepthType

_lh = get_logger(__name__)


def remark_fastq_single_end(
        src_fastq_file_path: str,
        dst_fastq_file_writer: FastqWriter,
        transcript_id: str,
        transcript_depth: float,
        simulator_name: str,
        truncate_ratio_3p: float,
        truncate_ratio_5p: float
) -> Tuple[int, int]:
    """
    Re-mark all seq_id in FASTQ files with truncation.

    - R: ``{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}``

    :return : Tuple of (number of actually processed reads, number of actually processed bases)

    .. versionadded:: 3.1.5
    """
    num_of_reads = 0
    num_of_bases = 0
    for fastq_record in FastqIterator(src_fastq_file_path, show_tqdm=False):
        sequence, quality = fastq_record.sequence, fastq_record.quality
        seq_len = len(sequence)
        truncate_len_3p = int(seq_len * truncate_ratio_3p)
        truncate_len_5p = int(seq_len * truncate_ratio_5p)
        sequence = sequence[truncate_len_5p:seq_len - truncate_len_3p]
        quality = quality[truncate_len_5p:seq_len - truncate_len_3p]
        new_record = FastqRecord(
            seq_id=f"{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}",
            sequence=sequence,
            quality=quality
        )
        dst_fastq_file_writer.write(new_record)
        num_of_reads += 1
        num_of_bases += len(quality)
    return num_of_reads, num_of_bases


def remark_fastq_pair_end(
        src_fastq_1_file_path: str,
        src_fastq_2_file_path: str,
        dst_fastq_1_file_writer: FastqWriter,
        dst_fastq_2_file_writer: FastqWriter,
        transcript_id: str,
        transcript_depth: float,
        simulator_name: str
) -> Tuple[int, int]:
    """
    Re-mark all seq_id in FASTQ files, return number of reads

    - R1: ``{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}/1``
    - R2: ``{transcript_id}:{num_of_reads}:{transcript_depth}:{simulator_name}/2``

    :return : Tuple of (number of actually processed reads, number of actually processed bases)

    .. versionadded:: 3.1.5
    """
    num_of_reads = 0
    num_of_bases = 0
    for fastq_record_1, fastq_record_2 in zip(
            FastqIterator(src_fastq_1_file_path, show_tqdm=False),
            FastqIterator(src_fastq_2_file_path, show_tqdm=False)
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
        dst_fastq_1_file_writer.write(new_fastq_record_1)
        dst_fastq_2_file_writer.write(new_fastq_record_2)
        num_of_reads += 1
        num_of_bases += len(new_fastq_record_1)
        num_of_bases += len(new_fastq_record_2)
    return num_of_reads, num_of_bases


class AssembError(ValueError):
    """
    TODO docs
    
    .. versionadded:: 3.1.5
    """
    ...


class AssemblerType(ABC, threading.Thread):
    """
    The Assembler that assembles sequenced transcripts of each isoform into one (SE) or two (PE) files.

    .. versionadded:: 3.1.5
    """

    def __init__(self, **kwargs) -> None:
        super().__init__()

    @abstractmethod
    def run(self) -> None:
        """Execute the assembler."""
        raise NotImplementedError

    @abstractmethod
    def add_transcript_id(self, transcript_id: str) -> None:
        """
        Add transcript_id to assembler.

        :param transcript_id: Transcript ID of an finished LLRG run.
        """
        raise NotImplementedError

    @abstractmethod
    def terminate(self) -> None:
        """
        Mark that no further Transcript ID would be added.

        Assembler would terminate after processing all added Transcript IDs.
        """
        raise NotImplementedError


class BaseAssembler(AssemblerType, ABC):
    """
    TODO docs
    
    .. versionadded:: 3.1.5
    """
    _transcript_ids_pending: List[str]
    _should_stop: bool
    _depth: DepthType
    _output_fastq_prefix: str
    _simulator_name: str
    _input_transcriptome_fasta_dir: str
    _input_fastq_dir: str

    def __init__(
            self,
            *,
            depth_data: DepthType,
            output_fastq_prefix: str,
            simulator_name: str,
            input_transcriptome_fasta_dir: str,
            input_fastq_dir: Optional[str] = None,
            **kwargs
    ):
        """
        Assembler Initialization.

        :param depth_data: Isoform-Level Depth Data.
        :param output_fastq_prefix: Prefix of output FASTQ used in LLRG orchestrators.
        :param simulator_name: Name of simulator. Used in FASTQ SeqID.
        :param input_transcriptome_fasta_dir: Directory of Transcriptome FASTA.
        :param kwargs: Other Miscellaneous arguments for compatibility. See subclass for details.
        """
        super().__init__()
        if input_fastq_dir is None:
            input_fastq_dir = output_fastq_prefix + ".d"
        self._transcript_ids_pending = []
        self._should_stop = False
        self._depth = depth_data
        self._output_fastq_prefix = output_fastq_prefix
        self._simulator_name = simulator_name
        self._input_transcriptome_fasta_dir = input_transcriptome_fasta_dir
        self._input_fastq_dir = input_fastq_dir

    def add_transcript_id(self, transcript_id: str):
        self._transcript_ids_pending.append(transcript_id)
        _lh.debug("ASSEMB: %s QUEUED", transcript_id)

    def terminate(self):
        _lh.info("ASSEMB: Termination signal received")
        self._should_stop = True

    @property
    def n_pending(self) -> int:
        return len(self._transcript_ids_pending)


class AssembleSingleEnd(BaseAssembler):
    """
    Assembler for single-end reads.
    
    .. versionadded:: 3.1.5
    """
    _truncate_ratio_3p: float
    _truncate_ratio_5p: float

    def __init__(
            self,
            *,
            depth_data: DepthType,
            output_fastq_prefix: str,
            simulator_name: str,
            input_transcriptome_fasta_dir: str,
            input_fastq_dir: Optional[str] = None,
            truncate_ratio_3p: float = 0,
            truncate_ratio_5p: float = 0,
            **kwargs
    ):
        """
        Assembler Initialization.

        :param depth_data: Isoform-Level Depth Data.
        :param output_fastq_prefix: Prefix of output FASTQ used in LLRG orchestrators.
        :param simulator_name: Name of simulator. Used in FASTQ SeqID.
        :param input_transcriptome_fasta_dir: Directory of Transcriptome FASTA.
        :param truncate_ratio_3p: Ratio of 3 prime truncation, range from 0 (no truncation) to 1 (truncate all).
        :param truncate_ratio_5p: Ratio of 5 prime truncation, range from 0 (no truncation) to 1 (truncate all).
        """
        super().__init__(
            depth_data=depth_data,
            output_fastq_prefix=output_fastq_prefix,
            simulator_name=simulator_name,
            input_transcriptome_fasta_dir=input_transcriptome_fasta_dir,
            input_fastq_dir=input_fastq_dir,
            **kwargs
        )
        self._truncate_ratio_3p = truncate_ratio_3p
        self._truncate_ratio_5p = truncate_ratio_5p

    def run(self):
        _lh.info("ASSEMB: START")
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

            def _assemb(transcript_id: str):
                _lh.debug("ASSEMB: %s START", transcript_id)
                this_fasta_path = os.path.join(self._input_transcriptome_fasta_dir, transcript_id + ".fa")
                this_fastq_basename = os.path.join(self._input_fastq_dir, transcript_id)
                this_fastq_path = this_fastq_basename + ".fq"
                if not file_system.file_exists(this_fastq_path):
                    _lh.warning("ASSEMB: %s: Skipped non-expressing transcript (FASTQ not exist)", transcript_id)
                    raise AssembError
                if not file_system.file_exists(this_fasta_path):
                    _lh.warning("ASSEMB: %s: Skipped none-existing transcript (FASTA not exist)", transcript_id)
                    raise AssembError
                try:
                    transcribed_length = FastaViewFactory(
                        filename=this_fasta_path,
                        read_into_memory=True,
                        show_tqdm=False
                    ).get_chr_length(transcript_id)
                except KeyError:
                    _lh.warning("ASSEMB: %s: Skipped none-existing transcript (FASTA Parsing Error)", transcript_id)
                    raise AssembError
                transcript_depth = self._depth[transcript_id]
                try:
                    num_of_reads, num_of_bases = remark_fastq_single_end(
                        src_fastq_file_path=this_fastq_path,
                        dst_fastq_file_writer=writer,
                        transcript_id=transcript_id,
                        transcript_depth=transcript_depth,
                        simulator_name=self._simulator_name,
                        truncate_ratio_3p=self._truncate_ratio_3p,
                        truncate_ratio_5p=self._truncate_ratio_5p
                    )
                except MisFormattedFastqRecordError:
                    _lh.warning("ASSEMB: %s: Skipped none-existing transcript (FASTQ Parsing Error)", transcript_id)
                    return
                stats_writer.write("\t".join((
                    transcript_id,  # "TRANSCRIPT_ID",
                    str(transcript_depth),  # "INPUT_DEPTH",
                    str(num_of_reads),  # "SIMULATED_N_OF_READS",
                    str(num_of_bases),  # "SIMULATED_N_OF_BASES",
                    str(transcribed_length),  # "TRANSCRIBED_LENGTH",
                    str(num_of_bases / transcribed_length)  # "SIMULATED_DEPTH"
                )) + "\n")
                _lh.debug("ASSEMB: %s FIN", transcript_id)
                stats_writer.flush()

            while not self._should_stop:
                if len(self._transcript_ids_pending) > 0:
                    try:
                        _assemb(self._transcript_ids_pending.pop())
                    except AssembError:
                        pass
                else:
                    time.sleep(0.01)
            _lh.info("ASSEMB: Post TERM -- Remaining %d transcripts", len(self._transcript_ids_pending))
            while len(self._transcript_ids_pending) > 0:
                try:
                    _assemb(self._transcript_ids_pending.pop())
                except AssembError:
                    pass
        _lh.info("ASSEMB: FIN")


class AssemblePairEnd(BaseAssembler):
    """
    Assembler for pair-end reads.
    
    .. versionadded:: 3.1.5
    """

    def run(self):
        def _assemb(transcript_id: str):
            _lh.debug("ASSEMB: %s START", transcript_id)
            this_fasta_path = os.path.join(self._input_transcriptome_fasta_dir, transcript_id + ".fa")
            this_fastq_basename = os.path.join(output_fastq_dir, transcript_id)
            this_fastq_r1_path = this_fastq_basename + "_1.fq"
            this_fastq_r2_path = this_fastq_basename + "_2.fq"
            if not file_system.file_exists(this_fastq_r1_path):
                _lh.warning("ASSEMB: %s: Skipped non-expressing transcript (FASTQ 1 not exist)", transcript_id)
                raise AssembError
            if not file_system.file_exists(this_fastq_r2_path):
                _lh.warning("ASSEMB: %s: Skipped non-expressing transcript (FASTQ 2 not exist)", transcript_id)
                raise AssembError
            if not file_system.file_exists(this_fasta_path):
                _lh.warning("ASSEMB: %s: Skipped none-existing transcript (FASTA not exist)", transcript_id)
                raise AssembError
            try:
                transcribed_length = FastaViewFactory(
                    filename=this_fasta_path,
                    read_into_memory=True,
                    show_tqdm=False
                ).get_chr_length(transcript_id)
            except KeyError:
                _lh.warning("ASSEMB: %s: Skipped none-existing transcript (FASTA Parsing Error)", transcript_id)
                raise AssembError
            transcript_depth = self._depth[transcript_id]
            try:
                num_of_reads, num_of_bases = remark_fastq_pair_end(
                    src_fastq_1_file_path=this_fastq_r1_path,
                    src_fastq_2_file_path=this_fastq_r2_path,
                    dst_fastq_1_file_writer=writer1,
                    dst_fastq_2_file_writer=writer2,
                    transcript_id=transcript_id,
                    transcript_depth=transcript_depth,
                    simulator_name=self._simulator_name
                )
            except MisFormattedFastqRecordError:
                _lh.warning("ASSEMB: %s: Skipped none-existing transcript (FASTQ Parsing Error)", transcript_id)
                return
            stats_writer.write("\t".join((
                transcript_id,  # "TRANSCRIPT_ID",
                str(transcript_depth),  # "INPUT_DEPTH",
                str(num_of_reads),  # "SIMULATED_N_OF_READS",
                str(num_of_bases),  # "SIMULATED_N_OF_BASES",
                str(transcribed_length),  # "TRANSCRIBED_LENGTH",
                str(num_of_bases / transcribed_length)  # "SIMULATED_DEPTH"
            )) + "\n")
            _lh.debug("ASSEMB: %s FIN", transcript_id)
            stats_writer.flush()

        _lh.info("ASSEMB: START")
        output_fastq_dir = self._output_fastq_prefix + ".d"
        with FastqWriter(self._output_fastq_prefix + "_1.fq") as writer1, \
                FastqWriter(self._output_fastq_prefix + "_2.fq") as writer2, \
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
                if len(self._transcript_ids_pending) > 0:
                    try:
                        _assemb(self._transcript_ids_pending.pop())
                    except AssembError:
                        pass
                else:
                    time.sleep(0.01)
            _lh.info("ASSEMB: Post TERM -- Remaining %d transcripts", len(self._transcript_ids_pending))
            while len(self._transcript_ids_pending) > 0:
                try:
                    _assemb(self._transcript_ids_pending.pop())
                except AssembError:
                    pass

        _lh.info("ASSEMB: FIN")


class AssembleDumb(BaseAssembler):
    """
    Assembler that does not perform anything.
    
    .. versionadded:: 3.1.5
    """

    def run(self):
        _lh.info("ASSEMB: START")
        while not self._should_stop:
            time.sleep(0.01)
        _lh.info("ASSEMB: FIN")


def generate_callback(
        assembler: AssemblerType,
        transcript_id: str
) -> Callable[[threading.Thread], None]:
    """
    Generate callback function for assemblers.

    This function is used by schedulers.
    It would invoke :py:func:`AssemblerType.add_transcript_id` when an instance of LLRG is finished.

    :param assembler: Targeted assembler.
    :param transcript_id: Source Transcript ID.
    :return: Generated callback.
    
    .. versionadded:: 3.1.5
    """
    return lambda _: assembler.add_transcript_id(transcript_id)


def patch_frontend_parser_public(
        parser: argparse.ArgumentParser,
        llrg_name: str,
        default_llrg_executable_name: Optional[str] = None
) -> argparse.ArgumentParser:
    """
    Patch argument parser with commonly-used options.

    :param llrg_name: Name of LLRG.
    :param default_llrg_executable_name: Default LLRG Executable Name. None if is function-based.
    :param parser: Source parser.
    :return: Patched parser.
    
    .. versionadded:: 3.1.5
    """
    parser.add_argument(
        '-F',
        '--fastas',
        required=True,
        help="Directory of transcribed cDNA sequences in FASTA format from `transcribe` step",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-j',
        '--jobs',
        required=False,
        help="Number of LLRGs to be executed in parallel",
        nargs='?',
        type=int,
        action='store',
        default=multiprocessing.cpu_count()
    )
    parser.add_argument(
        '--simulator_name',
        required=False,
        help="Custom simulator name. Used in FASTQ tags\n"
             "This step is done in assemble, so if --not_perform_assemble is set, this option would be useless.",
        nargs='?',
        type=str,
        action='store',
        default=None
    )
    if default_llrg_executable_name is not None:
        parser.add_argument(
            '-e',
            '--llrg_executable_path',
            required=False,
            help=f"Executable name or absolute path of {llrg_name}",
            nargs='?',
            type=str,
            action='store',
            default=default_llrg_executable_name
        )
    return parser


def patch_frontend_parser_sc_rna_seq(
        parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Patch argument parser for Single-Cell RNA-Seq

    :param parser: Source parser.
    :return: Patched parser.
    
    .. versionadded:: 3.1.5
    """
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        help="Path to Isoform-Level Depth directory generated by `generate_barcoded_isoform_replicates`",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Path to output Directory.",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def patch_frontend_parser_bulk_rna_seq_assemble(
        parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    TODO docs
    
    .. versionadded:: 3.1.5
    """
    return parser


def patch_frontend_parser_bulk_rna_seq(
        parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Patch argument parser for Bulk-Cell RNA-Seq

    :param parser: Source parser.
    :return: Patched parser.
    
    .. versionadded:: 3.1.5
    """
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        help="Path to input Isoform-Level Depth TSV generated by `generate_depth_v2` or `generate_isoform_depth` step",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output transcript prefix. "
             "The output file would be {out}.fq for single-end "
             "and {out}_1.fq, {out}_2.fq for pair end.\n"
             "If --not_perform_assemble is set, would NOT generate FASTQ files but a {out}.d unassembled directory",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--not_perform_assemble',
        help="Do NOT assemble the output of each isoforms into one file.",
        action='store_true'
    )
    return parser


def patch_frontend_parser_tgs(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Patch argument parser with public option and
    ``--truncate_ratio_3p`` and ``--truncate_ratio_5p`` option.

    :param parser: Source parser.
    :return: Patched parser.
    
    .. versionadded:: 3.1.5
    """
    parser.add_argument(
        '--truncate_ratio_3p',
        required=False,
        help="[Single End TGS Only] Ratio of 3 prime truncation, range from 0 (no truncation) to 1 (truncate all).\n"
             "This step is done in assemble, so if --not_perform_assemble is set, this option would be useless.",
        nargs='?',
        type=float,
        action='store',
        default=0.0
    )
    parser.add_argument(
        '--truncate_ratio_5p',
        required=False,
        help="[Single End TGS Only] Ratio of 5 prime truncation, range from 0 (no truncation) to 1 (truncate all).\n"
             "This step is done in assemble, so if --not_perform_assemble is set, this option would be useless.",
        nargs='?',
        type=float,
        action='store',
        default=0.0
    )
    return parser


def enhanced_which(path_or_filename: str) -> str:
    """
    Enhanced :py:func:`shutil.which` method that would raise error on failure.

    :param path_or_filename: Path or filename you wish to resolve.
    :return: Absolute path found by :py:func:`shutil.which`.
    :raise FileNotFoundError: On failure.
    
    .. versionadded:: 3.1.5
    """
    resolved_path = shutil.which(path_or_filename)
    if resolved_path is None:
        raise FileNotFoundError(f"File {path_or_filename} not found!")
    else:
        path_or_filename = resolved_path
    return path_or_filename
