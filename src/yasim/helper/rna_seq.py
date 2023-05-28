import glob
import os
import sys
import time
from collections import defaultdict

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio import file_system
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor
from labw_utils.typing_importer import Tuple, Mapping, Any, Type, Optional, Dict, List
from yasim.helper.depth_io import DepthType, read_depth, DepthParsingException
from yasim.helper.llrg import enhanced_which, AssembleDumb, AssemblePairEnd, AssembleSingleEnd, generate_callback
from yasim.llrg_adapter import LLRGInitializationException, BaseLLRGAdapter

_lh = get_logger(__name__)


def validate_adapter_args(
        adapter_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        llrg_executable_path: Optional[str] = None
) -> Mapping[str, Any]:
    _lh.info("Validating LLRG Adapter parameters...")
    # Validate LLRG Adapter Params
    if llrg_executable_path is not None:
        try:
            llrg_executable_path = enhanced_which(llrg_executable_path)
        except FileNotFoundError:
            _lh.error("LLRG not found at %s", llrg_executable_path)
            sys.exit(1)
    adapter_args = dict(adapter_args)
    try:
        adapter_args_update = adapter_class.validate_params(
            **adapter_args
        )
    except LLRGInitializationException as e:
        _lh.error("Exception %s caught at validating LLRG adapter arguments.", str(e))
        sys.exit(1)
    adapter_args.update(adapter_args_update)
    return adapter_args


def run_rna_seq(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        depth_data: DepthType,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        llrg_executable_path: Optional[str] = None,
        not_perform_assemble: bool = False,
        show_tqdm: bool = True
):
    _lh.info("Creating multiprocessing pool and assembler...")
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs,
        delete_after_finish=False,
        show_tqdm=show_tqdm
    )
    depth_info: List[Tuple[float, str, str]] = []
    """
    Depth information used by LLRG frontend interfaces.
    
    They are: [depth, transcript_id, filename]
    """

    for transcript_id, transcript_depth in depth_data.items():
        filename = os.path.join(transcriptome_fasta_dir, transcript_id + ".fa")
        if not file_system.file_exists(filename):
            _lh.warning("FASTA of Isoform %s not exist!", transcript_id)
            continue
        if transcript_depth <= 0:
            _lh.warning("Depth of Isoform %s too low!", transcript_id)
            continue
        depth_info.append((transcript_depth, transcript_id, filename))

    # Get assembler arguments
    assembler_full_args = {
        "depth_data": depth_data,
        "output_fastq_prefix": output_fastq_prefix,
        "simulator_name": simulator_name,
        "input_transcriptome_fasta_dir": transcriptome_fasta_dir,
    }
    assembler_full_args.update(assembler_args)

    # Create assembler
    if not_perform_assemble:
        assembler_class = AssembleDumb
    elif is_pair_end:
        assembler_class = AssemblePairEnd
    else:
        assembler_class = AssembleSingleEnd
    assembler = assembler_class(**assembler_full_args)
    assembler.start()

    _lh.info("Starting simulation and assembly...")
    for transcript_depth, transcript_id, transcript_filename in tqdm(
            iterable=depth_info,
            desc="Submitting jobs..."
    ):
        try:
            sim_thread = adapter_class(
                src_fasta_file_path=transcript_filename,
                dst_fastq_file_prefix=os.path.join(output_fastq_dir, transcript_id),
                depth=transcript_depth,
                llrg_executable_path=llrg_executable_path,  # TODO: Fix this.
                is_trusted=True,
                **adapter_args
            )
        except LLRGInitializationException as e:
            _lh.warning("Exception %s caught at initialization time", str(e))
            continue
        simulating_pool.append(sim_thread, callback=generate_callback(assembler, transcript_id))

    _lh.info("Jobs submitted, waiting...")
    simulating_pool.start()
    simulating_pool.join()

    _lh.info("Jobs finished, terminating assembler...")
    assembler.terminate()
    _lh.info("Assembler termination signal sent, waiting...")
    while assembler.is_alive():
        _lh.info("%s -- PENDING: %d", output_fastq_prefix,  assembler.n_pending)
        time.sleep(1.0)
    assembler.join()
    _lh.info("Assembler finished, retrieving error reports...")

    # Calculating number of exceptions
    exception_dict = defaultdict(lambda: 0)
    for job in simulating_pool.iter_finished_jobs():
        job: BaseLLRGAdapter
        exception_dict[job.exception] += 1
    return exception_dict


def bulk_rna_seq_frontend(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        depth_file_path: str,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        llrg_executable_path: Optional[str] = None,
        not_perform_assemble: bool = False
) -> int:
    """
    :return: Exit Value
    """
    adapter_args = validate_adapter_args(
        adapter_args=adapter_args,
        adapter_class=adapter_class,
        llrg_executable_path=llrg_executable_path
    )
    try:
        depth_data = read_depth(depth_file_path)
    except DepthParsingException:
        _lh.error(f"Failed to parse depth file {depth_file_path}")
        return 1

    exception_dict = run_rna_seq(
        transcriptome_fasta_dir=transcriptome_fasta_dir,
        output_fastq_prefix=output_fastq_prefix,
        depth_data=depth_data,
        jobs=jobs,
        simulator_name=simulator_name,
        adapter_args=adapter_args,
        assembler_args=assembler_args,
        adapter_class=adapter_class,
        is_pair_end=is_pair_end,
        llrg_executable_path=llrg_executable_path,
        not_perform_assemble=not_perform_assemble,
        show_tqdm=True
    )
    _lh.info(f"Status of errors: {dict(exception_dict)}")
    _lh.info("Simulation finished successfully")
    return 0


def sc_rna_seq_frontend(
        transcriptome_fasta_dir: str,
        output_fastq_dir_path: str,
        depth_dir_path: str,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        llrg_executable_path: Optional[str] = None,
        not_perform_assemble: bool = False
) -> int:
    """
    :return: Exit Value
    """
    adapter_args = validate_adapter_args(
        adapter_args=adapter_args,
        adapter_class=adapter_class,
        llrg_executable_path=llrg_executable_path
    )

    barcode_depth_data_dict: Dict[str, DepthType] = {}
    for depth_file_path in tqdm(glob.glob(os.path.join(depth_dir_path, "*")), desc="Validating depth input"):
        try:
            depth_data = read_depth(depth_file_path)
        except DepthParsingException:
            _lh.error(f"Failed to parse depth file {depth_file_path}")
            return 1
        barcode = os.path.basename(depth_file_path).split(".")[0]
        if barcode in barcode_depth_data_dict:
            _lh.error(f"Duplicated barcode {barcode}")
            return 1
        barcode_depth_data_dict[barcode] = depth_data

    full_exception_dict = defaultdict(lambda: 0)
    for barcode, depth_data in barcode_depth_data_dict.items():
        exception_dict = run_rna_seq(
            transcriptome_fasta_dir=transcriptome_fasta_dir,
            output_fastq_prefix=os.path.join(output_fastq_dir_path, barcode),
            depth_data=depth_data,
            jobs=jobs,
            simulator_name=simulator_name,
            adapter_args=adapter_args,
            assembler_args=assembler_args,
            adapter_class=adapter_class,
            is_pair_end=is_pair_end,
            llrg_executable_path=llrg_executable_path,
            not_perform_assemble=not_perform_assemble,
            show_tqdm=False
        )
        for k, v in exception_dict.items():
            full_exception_dict[k] += v

    _lh.info(f"Status of errors: {dict(full_exception_dict)}")
    _lh.info("Simulation finished successfully")
    return 0
