import os
from collections import defaultdict
from typing import Mapping, Any, Type

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor
from yasim.helper.depth import DepthType
from yasim.helper.llrg import AssembleDumb, pair_depth_info_with_transcriptome_fasta_filename, \
    generate_callback, AssembleSingleEnd, AssemblePairEnd
from yasim.llrg_adapter import BaseLLRGAdapter, LLRGInitializationException

_lh = get_logger(__name__)


def abstract_simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        depth: DepthType,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        not_perform_assemble: bool = False
):
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs,
        delete_after_finish=False
    )
    depth_info = list(pair_depth_info_with_transcriptome_fasta_filename(transcriptome_fasta_dir, depth))

    # Get assembler arguments
    assembler_full_args = {
        "depth_data": depth,
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

    # Start simulation and assembly.
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        try:
            sim_thread = adapter_class(
                src_fasta_file_path=transcript_filename,
                dst_fastq_file_prefix=os.path.join(output_fastq_dir, transcript_id),
                depth=transcript_depth,
                is_trusted=True,
                **adapter_args
            )
        except LLRGInitializationException as e:
            _lh.error("Exception %s caught at initialization time", str(e))
            continue
        simulating_pool.append(sim_thread, callback=generate_callback(assembler, transcript_id))
    simulating_pool.start()
    simulating_pool.join()

    # Terminate assemblers.
    assembler.terminate()
    assembler.join()

    # Calculating number of exceptions
    exception_dict = defaultdict(lambda: 0)
    for job in simulating_pool.iter_finished_jobs():
        job: BaseLLRGAdapter
        exception_dict[job.exception] += 1
    _lh.info(f"Status of errors: {dict(exception_dict)}")
