import os
from collections import defaultdict
from typing import Mapping, Any, Type

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor
from yasim.helper.depth import read_depth
from yasim.helper.llrg import AssembleDumb, pair_depth_info_with_transcriptome_fasta_filename, \
    generate_callback, AssembleSingleEnd, AssemblePairEnd, enhanced_which
from yasim.llrg_adapter import BaseLLRGAdapter, LLRGInitializationException

_lh = get_logger(__name__)


def abstract_simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        depth_file_path: str,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        llrg_executable_path: str,
        not_perform_assemble: bool
) -> int:
    """

    :param transcriptome_fasta_dir:
    :param output_fastq_prefix:
    :param depth_file_path:
    :param jobs:
    :param simulator_name:
    :param adapter_args:
    :param assembler_args:
    :param adapter_class:
    :param is_pair_end:
    :param llrg_executable_path:
    :param not_perform_assemble:
    :return: Exit Value
    """
    _lh.info("Validating LLRG Adapter parameters...")
    # Validate LLRG Adapter Params
    adapter_args = dict(adapter_args)
    try:
        adapter_args_update = adapter_class.validate_params(
            **adapter_args
        )
    except LLRGInitializationException as e:
        _lh.error("Exception %s caught at validating LLRG adapter arguments.", str(e))
        return 1
    adapter_args.update(adapter_args_update)
    try:
        llrg_executable_path = enhanced_which(llrg_executable_path)
    except FileNotFoundError:
        _lh.error("LLRG not found at %s", llrg_executable_path)
        return 1

    depth_data = read_depth(depth_file_path)

    _lh.info("Creating multiprocessing pool and assembler...")
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs,
        delete_after_finish=False
    )
    depth_info = list(pair_depth_info_with_transcriptome_fasta_filename(transcriptome_fasta_dir, depth_data))

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
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        try:
            sim_thread = adapter_class(
                src_fasta_file_path=transcript_filename,
                dst_fastq_file_prefix=os.path.join(output_fastq_dir, transcript_id),
                depth=transcript_depth,
                llrg_executable_path=llrg_executable_path,
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

    _lh.info("Jobs finished, cleaning up...")
    assembler.terminate()
    assembler.join()

    # Calculating number of exceptions
    exception_dict = defaultdict(lambda: 0)
    for job in simulating_pool.iter_finished_jobs():
        job: BaseLLRGAdapter
        exception_dict[job.exception] += 1
    _lh.info(f"Status of errors: {dict(exception_dict)}")
    _lh.info("Simulation finished successfully")
    return 0
